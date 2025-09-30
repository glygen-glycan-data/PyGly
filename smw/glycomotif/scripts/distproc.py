#!/bin/env python3.12

import random
import time
import sys
import socket
import queue
import secrets
import os
import os.path
import hashlib
import math
import multiprocessing
import traceback
import threading
import subprocess
import signal
from collections import defaultdict, deque
from tqdm import tqdm

from multiprocessing.managers import SyncManager

class DistributedProcessing(object):
    def __init__(self,target=None,host=None,port=None,secret=None,verbose=False):
        self.hostname = socket.gethostname()
        if target:
            self.target = target
        self.verbose = verbose
        self.serialproc = False
        self.manager = None
        self.worker_procs = []
        self.procs = []
        self.incleanup = False
        if host:
            self.host = host
        else:
            self.host = self.hostname
        if port:
            self.port = int(port)
        else:
            self.port = self.random_port()
        if secret:
            self.secret = secret.encode()
        else:
            self.secret = self.random_secret()

    @staticmethod
    def random_port(seed_string=None):
        if seed_string:
            return (60900 + int(hashlib.sha256(seed_string.encode()).hexdigest(),16)%100)
        return random.randint(60900,60999)

    @staticmethod
    def random_secret(seed_string=None,length=16):
        if seed_string:
            return hashlib.sha256(seed_string.encode()).hexdigest()[:length].encode()
        return ("".join([random.choice('0123456789abcdef') for _ in range(length)])).encode()

    def make_server_manager(self):

        class JobQueueManager(SyncManager):
            pass

        self._tasks =  queue.Queue()
        self._results = queue.Queue()
        self._workers = queue.Queue()
        self._manager = queue.Queue()
        self._data = dict()

        JobQueueManager.register('get_task_queue', callable=lambda: self._tasks)
        JobQueueManager.register('get_result_queue', callable=lambda: self._results)
        JobQueueManager.register('get_worker_queue', callable=lambda: self._workers)
        JobQueueManager.register('get_manager_queue', callable=lambda: self._manager)
        JobQueueManager.register('get_shared_data', callable=lambda: self._data)
        
        self.register_cleanup()
        self.manager = JobQueueManager(address=("", self.port), authkey=self.secret)
        self.manager.start()
        if self.verbose == True:
            print('Server started at port %s (secret: %s)' % (self.port, self.secret.decode()), file=sys.stderr)
        self.tasks = self.manager.get_task_queue()
        self.results = self.manager.get_result_queue()
        self.worker_messages = self.manager.get_worker_queue()
        self.manager_messages = self.manager.get_manager_queue()
        self.shared_data = self.manager.get_shared_data()
        self.njobs = 0

    def make_client_manager(self):
        class ServerQueueManager(SyncManager):
            pass

        ServerQueueManager.register('get_task_queue')
        ServerQueueManager.register('get_result_queue')
        ServerQueueManager.register('get_worker_queue')
        ServerQueueManager.register('get_manager_queue')
        ServerQueueManager.register('get_shared_data')

        if self.verbose == True:
            print('Client attempting connection to %s:%s (%s)' % (self.host, self.port, self.secret.decode()), file=sys.stderr)
        self.register_cleanup()
        self.manager = ServerQueueManager(address=(self.host, self.port), authkey=self.secret)
        ntries = 4
        for i in range(ntries):
            try:
                self.manager.connect()
                break
            except IOError:
                if i == (ntries - 1):
                    raise
                print('Client failed to connect, attempt %d'%(i+1,),file=sys.stderr)
                time.sleep(5)
        if self.verbose == True:
            print('Client connected to %s:%s (%s)' % (self.host, self.port, self.secret.decode()), file=sys.stderr)
        self.tasks = self.manager.get_task_queue()
        self.results = self.manager.get_result_queue()
        self.worker_messages = self.manager.get_worker_queue()
        self.manager_messages = self.manager.get_manager_queue()
        self.shared_data = self.manager.get_shared_data()

    def start_workers(self,ncpus):
        self.worker_procs = []
        pid = os.getpid()
        for i in range(ncpus):
            workerid = "%s:%s"%(pid,i+1)
            proc = multiprocessing.Process(target=self.worker,args=(workerid,),daemon=True)
            self.worker_procs.append(proc)
            proc.start()

    def put_task(self,task_index,task):
        try:
            self.tasks.put(dict(task_index=task_index,task=task))
        except (BrokenPipeError,EOFError):
            pass

    def get_task(self):
        try:
            return self.tasks.get()
        except (BrokenPipeError,EOFError):
            pass
        return None

    def put_result(self,worker_index,task,task_index,elapsed,result):
        try:
            self.results.put(dict(status="RESULT",hostname=self.hostname,worker_index=worker_index,
                                  task=task,task_index=task_index,runtime=elapsed,result=result,
                                  stage=task.get('__stage__')))
        except (BrokenPipeError,EOFError):
            pass

    def put_error(self,worker_index,task,task_index,elapsed,excep):
        try:
            self.results.put(dict(status="ERROR",hostname=self.hostname,worker_index=worker_index,
                                  task=task,task_index=task_index,runtime=elapsed,
                                  traceback=traceback.format_exception(*excep),
                                  stage=task.get('__stage__')))
        except (BrokenPipeError,EOFError):
            pass

    def get_result(self,timeout=-1):
        try:
            return self.results.get(timeout=timeout)
        except (queue.Empty,BrokenPipeError,EOFError):
            return None

    def do_task(self,task,**kwargs):
        if self.target is None:
            raise NotImplemented("Neither target nor derived class do_task method defined.")
        if task.get('__stage__') is not None:
            return self.target[task['__stage__']](task['__task__'],**kwargs)
        return self.target(task['__task__'],**kwargs)

    def init(self):
        return

    def heartbeat(self,worker_index):
        while True:
            self.worker_messages.put(("HEARTBEAT",worker_index))
            time.sleep(60)

    def worker(self,worker_index):
        self.worker_messages.put(("WORKERID",worker_index))
        t=threading.Thread(target=self.heartbeat,args=(worker_index,))
        t.daemon = True
        t.start()
        init_called = False
        while True:
            task = self.get_task()
            if task is None:
                break
            task_index = task['task_index']
            task = task['task']
            if task is None:
                break
            self.worker_messages.put(("TASKID",task_index,worker_index))
            if not init_called:
                self.init()
                init_called = True
            try:
                start = time.time()
                result = self.do_task(task,hostname=self.hostname,worker_index=worker_index,task_index=task_index,shared_data=self.shared_data)
            except Exception:
                self.put_error(worker_index,task,task_index,int(round(time.time()-start,0)),sys.exc_info())
            else:
                self.put_result(worker_index,task,task_index,int(round(time.time()-start,0)),result)
        return

    def wait_workers(self):
        for p in self.worker_procs:
            p.join()

    def shutdown(self):
        self.manager.shutdown()

    def client(self,ncpus):
        assert ncpus > 0
        self.make_client_manager()
        self.start_workers(ncpus)
        self.wait_workers()

    def cleanup(self,*args):
        self.incleanup = True
        for p in self.worker_procs:
            while p.is_alive():
                p.kill()
        if self.manager is not None:
            self.manager.shutdown()
        os._exit(130)

    def register_cleanup(self):
        signal.signal(signal.SIGINT, self.cleanup)

    def server(self):
        self.make_server_manager()
        return self

    def procspec(self,arg):
        worker_args = []
        for argi in sys.argv[1:]:
            if argi == arg:
                worker_args.append("__%(ncpus)s:%(server)s:%(port)s:%(secret)s__")
            else:
                worker_args.append(argi)
        procspec = {None: 0}
        for ps in arg.split(','):
            sps = ps.split(':')
            if len(sps) == 1:
                procspec[None] = int(sps[0])
            else:
                procspec[sps[0].strip()] = [ int(i) for i in sps[1:] ]
        for k,v in procspec.items():
            if not k:
                continue
            if k == "slurm":
                self.start_slurm_workers(v,worker_args)
            else:
                self.start_remote_workers(k,v,worker_args)
        return procspec[None]

    def start_remote_workers(self,worker,spec,worker_args):
        cmd = 'ssh -n -f %s nohup sh -c \\\'"cd %s; %s %s'%(worker,os.getcwd(),sys.executable,sys.argv[0])
        for arg in worker_args:
            cmd += " "+arg%dict(server=self.hostname,ncpus=spec[0],port=self.port,secret=self.secret.decode())
        cmd += " >/dev/null 2>&1 & \"\\\'"
        # print(cmd,file=sys.stderr)
        # p = subprocess.run(cmd,shell=True,check=True,stdin=subprocess.DEVNULL)
        p = subprocess.run(cmd,shell=True,check=True,stdin=subprocess.DEVNULL,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
        self.procs.append(p)

    def start_slurm_workers(self,spec,worker_args):
        if len(spec) >= 2:
            njobs = spec[0]
            ncpus = spec[1]
        else:
            njobs = spec[0]
            ncpus = 1
        sbatch = 'sbatch --cpus-per-task %s --output=/dev/null --array=1-%d'%(ncpus,njobs)
        cmd = '%s %s'%(sys.executable,os.path.abspath(sys.argv[0]))
        for arg in worker_args:
            cmd += " "+arg%dict(server=self.hostname,ncpus=ncpus,port=self.port,secret=self.secret.decode())
        stdinstr = "\n".join(map(str.lstrip,filter(None,"""
        #!/bin/sh
        srun %s
        """.splitlines())))
        # print(sbatch)
        # print(stdinstr%(cmd,))
        subprocess.run(sbatch,input=stdinstr%(cmd,),capture_output=True,text=True,check=True,shell=True)

    def startup(self,workers="",**shared_data):
        self.shared_data.update(shared_data)
        self.start_workers(self.procspec(workers))
        self.workerids = set()
        self.heartbeat = defaultdict(float)
        return self

    def execute(self,tasks,noshutdown=False,stage=None):
        self.alltasks = [ {'__stage__': stage, '__task__': t} for t in tasks ]
        self.noshutdown=noshutdown
        self.stage = stage        
        return self

    def ntasks(self):
        return len(self.alltasks)

    def tasksempty(self):
        try:
            return self.tasks.empty()
        except (BrokenPipeError,EOFError):
            pass
        return True

    def worker_messages_empty(self):
        try:
            return self.worker_messages.empty()
        except (BrokenPipeError,EOFError,AttributeError):
            pass
        return True

    def hms(self,seconds,lb=True):
        if lb:
            s = int(math.floor(seconds))
        else:
            s = int(math.ceil(seconds))
        h = int(math.floor(s/3600))
        s -= h*3600
        m = int(math.floor(s/60))
        s -= m*60
        return "%d:%02d:%02d"%(h,m,s)

    def hm(self,seconds,lb=True):
        if lb:
           m = int(math.floor(seconds/60))
        else:
           m = int(math.ceil(seconds/60))
        h = int(math.floor(m/60))
        m -= h*60
        return "%d:%02d"%(h,m)

    def update_progress(self,result):
        start = self.starttime
        now = time.time()
        done = len(self.donetasks)
        alltasks = len(self.alltasks)
        remain = alltasks-done
 
        result['done'] = done
        result['percent_done'] = 100*done/alltasks
        result['progress'] = "%s/%s (%.2f%%)"%(done, alltasks, 100*done/alltasks)
        if start is not None:
            elapsed = now-start
            result['elapsed'] = self.hms(elapsed,lb=True)
            result['remaining'] = self.hm(remain*elapsed/done,lb=False)
            result['taskspermin'] = "%.2f"%(60*done/elapsed,)
        else:
            result['elapsed'] = ""
            result['remaining'] = ""
            result['taskspermin'] = ""

    def __iter__(self):
        return self.iterresults()

    def iterresults(self):

        if self.verbose == True:
            print("Begin execution: %s tasks to complete."%(len(self.alltasks),),file=sys.stderr)

        for i,task in enumerate(self.alltasks):
            self.put_task(i+1,task)

        if self.verbose == True:
            print("%s tasks added to the task queue."%(len(self.alltasks),),file=sys.stderr)

        self.starttime = None
        self.donetasks = set()
        self.taskattempts = defaultdict(int)
        self.failedtasks = set()
        self.task2worker = dict()
        self.tasksemptytime = None
        while not self.tasksempty() or (len(self.donetasks) + len(self.failedtasks)) < len(self.alltasks):
 
            while not self.worker_messages_empty():
                msg = self.worker_messages.get()
                # print(msg,file=sys.stderr)
                if msg[0] == "WORKERID":
                    self.workerids.add(msg[1])
                elif msg[0] == "TASKID":
                    self.task2worker[msg[1]] = msg[2]
                    if self.starttime is None:
                        self.starttime = time.time()
                elif msg[0] == "HEARTBEAT":
                    self.heartbeat[msg[1]] = time.time()

            if self.tasksempty() and self.tasksemptytime is None:
                self.tasksemptytime = time.time()
                
            if self.tasksempty() and not self.incleanup and self.tasksemptytime is not None and (time.time() - self.tasksemptytime) > 30:
                for i,task in enumerate(self.alltasks):
                    taskid = (i+1)
                    if taskid not in self.donetasks and taskid not in self.failedtasks:
                        if taskid not in self.task2worker or (time.time() - self.heartbeat[self.task2worker[taskid]]) > 120:
                            print("Reqeueing missing task %d..."%(taskid,),file=sys.stderr)
                            self.put_task(taskid,task)

            result = self.get_result(15)
            if result is None:
                continue

            status = result['status']
            if status == "ERROR":
                print("Worker %s:%s: Task %s error...\n%s"%(result.get('hostname'),result.get('worker_index'),
                                                            result.get('task_index'),"".join(result.get('traceback',[]))),
                                                            end="",file=sys.stderr)
                taskid = result.get('task_index')
                self.taskattempts[taskid] += 1
                if self.taskattempts[taskid] < 3:
                    print("Reqeueing on error task %d..."%(taskid,),file=sys.stderr)
                    self.put_task(taskid,self.alltasks[taskid-1])
                else:
                    print("Failed task %d..."%(taskid,),file=sys.stderr)
                    self.failedtasks.add(taskid)
            elif status == "RESULT":
                taskid = result.get('task_index')
                if taskid not in self.donetasks and (not self.stage or self.stage == result.get('stage')):
                    self.donetasks.add(taskid)
                    self.update_progress(result)
                    yield result
            else:
                raise RuntimeError("Bad result status")

        if self.verbose == True:
            print("Task summary: %s tasks completed, %s tasks failed."%(len(self.donetasks),len(self.failedtasks)),file=sys.stderr)

        if not self.noshutdown:
            self.allshutdown()

    def allshutdown(self):
        while not self.worker_messages_empty():
            msg = self.worker_messages.get()
            if msg[0] == "WORKERID":
                self.workerids.add(msg[1])

        for i in range(len(self.workerids)):
            self.put_task(-1,None)
        
        self.wait_workers()  

        time.sleep(5)
        self.shutdown()

    def serial(self,**shared_data):
        self.shared_data = shared_data
        self.iterresults = self.serialiterresults
        self.serialproc = True
        return self

    def serialiterresults(self):
        pid = os.getpid()
        workerid = "%s"%(pid,)
        self.init()
        self.starttime = time.time()
        self.donetasks = set()
        for i,task in enumerate(self.alltasks):
            start = time.time()
            result = self.do_task(task,hostname=self.hostname,task_index=(i+1),worker_index=workerid,shared_data=self.shared_data)
            result = dict(status='RESULT',hostname=self.hostname,worker_index=workerid,task=task,task_index=(i+1),runtime=int(round(time.time()-start,0)),result=result)
            self.donetasks.add(i+1)
            self.update_progress(result)
            yield result           

    @staticmethod
    def add_arguments(parser,argname="workers"):
        helpstr = " ".join("""
                  Enable distributed processing: WORKERS = 
                  [<n0>,][remote1:<n1>[,remote2:<n2>,...]]. n0 is cpus on host node
                  (optional), ni is cpus on optional remote nodes. Remote
                  nodes yion, bion, proton, and maldi with 8 cpus each are
                  available. Use remote node name "slurm" for shared
                  computing resources, with 32 cpus available. Default:
                  serial, single cpu, processing.
        """.split())
        parser.add_argument('--'+argname, dest="__distproc__", metavar="WORKERS", type=str, default  = "", help = helpstr)
        return parser

    @staticmethod
    def parse_args(parser):
        args = parser.parse_args()
        distproc = args.__distproc__.strip()
        if distproc.startswith('__') and distproc.endswith('__'):
            return ('worker',distproc.strip('_'))
        if distproc == "":
            return None
        return ('manager',distproc)

    logtempl = "worker_id: %(hostname)s:%(worker_index)s task_id: %(task_index)s runtime: %(runtime)s progress: %(progress)s remaining: %(remaining)s"

    @staticmethod
    def process(workers,target,tasks,verbose=None,logtempl=logtempl,**shared_data):
        if workers is None: 
            # serial processing
            for result in DistributedProcessing(target=target,verbose=(verbose==True)).serial(**shared_data).execute(tasks):
                if verbose in (None,'TQDM',True):
                    print(logtempl%result,file=sys.stderr)
                yield result['result']

        elif workers[0] == "manager":
            # manager/server/hostnode
            p = DistributedProcessing(target=target,verbose=(verbose==True)).server()
            p.startup(workers=workers[1],**shared_data)
            nextoutput = time.time()-10
            iterator = p.execute(tasks)
            if verbose in ('TQDM',None):
                iterator = tqdm(iterator,ascii=True,total=p.ntasks())
            for result in iterator:
                if (verbose == True) and (time.time() > nextoutput):
                    print(logtempl%result,file=sys.stderr)
                    nextoutput = time.time()+15
                yield result['result']

        elif workers[0] == "worker":
            # worker
            ncpus,server,port,secret = workers[1].split(':')
            DistributedProcessing(target=target,host=server,port=port,secret=secret).client(int(ncpus))
            sys.exit(0)

    @staticmethod
    def start_if_worker(workers,target):
        if workers is not None and workers[0] == "worker":
            # worker
            ncpus,server,port,secret = workers[1].split(':')
            DistributedProcessing(target=target,host=server,port=port,secret=secret).client(int(ncpus))
            sys.exit(0)

    @staticmethod
    def stage_process_init(workers,targets,verbose=None,**shared_data):
        if workers is None: 
            # serial processing
            return DistributedProcessing(target=targets,verbose=verbose).serial(**shared_data)
        elif workers[0] == "manager":
            # manager/server/hostnode
            return DistributedProcessing(target=targets,verbose=verbose).server().startup(workers[1],**shared_data)
        elif workers[0] == "worker":
            # worker
            ncpus,server,port,secret = workers[1].split(':')
            DistributedProcessing(target=targets,host=server,port=port,secret=secret).client(int(ncpus))
            sys.exit(0)
        
    def stage_process(self,stage,tasks):
        nextoutput = time.time()-10
        iterator = self.execute(tasks,noshutdown=True,stage=stage)
        if self.verbose in ('TQDM',None):
            iterator = tqdm(iterator,desc="Stage %s"%(stage,),ascii=True,total=self.ntasks())
        for result in iterator:
            if (self.verbose == True) and (time.time() > nextoutput):
                print(self.logtempl%result,file=sys.stderr)
                nextoutput = time.time() + 15
            yield result['result']

    def stage_process_finish(self):
        if not self.serialproc:
            self.allshutdown()
                                                                                                         
def do_task(task,**kwargs):
    # print("Worker %s:%s: Task %s delay %s starting..."%(kwargs.get('hostname'),kwargs.get('worker_index'),
    #                                                     kwargs.get('task_index'),task),file=sys.stderr)
    # print("Shared data:",kwargs.get('shared_data'))
    time.sleep(task)
    assert(random.random() < 0.95)
    # print("Worker %s:%s: Task %s delay %s done..."%(kwargs.get('hostname'),kwargs.get('worker_index'),
    #                                                 kwargs.get('task_index'),task),file=sys.stderr)
    return task

def process_result(result,**kwargs):
    print("Result:",result,kwargs,file=sys.stderr)

if __name__ == "__main__":
    tasks = [ random.randint(0,20) for i in range(100) ]

    if len(sys.argv) <= 1:

        for result in DistributedProcessing(target=do_task).serial(tasks):
            process_result(**result)

    elif sys.argv[1] == "manager":

        p = DistributedProcessing(target=do_task,workerargs=["worker","%(ncpus)s","%(server)s"]).server()
        for result in p.execute(tasks,workers=sys.argv[2]):
            process_result(**result)
                        
    elif sys.argv[1] == "worker":
                        
        ncpus = int(sys.argv[2])
        DistributedProcessing(target=do_task,host=sys.argv[3]).client(ncpus)
