
from tempfile import mkstemp
import os,sys,re,os.path,time
from subprocess import Popen, PIPE, STDOUT
import threading

class PopenTimeout(threading.Thread):
    def __init__(self, *args, **kwargs):
        threading.Thread.__init__(self)
        self.timeout = None
        if kwargs.get('timeout') != None:
            self.timeout = kwargs['timeout']
	    del kwargs['timeout']
	self.args = args
	self.kwargs = kwargs
        self.p = None
	self.retval = None
        self.start()
        while self.p == None:
            time.sleep(.1)

    def run(self):
        self.p = Popen(*self.args,**self.kwargs)
        self.retval = self.p.wait()

    def wait(self):
        self.join(self.timeout)
        if self.is_alive():
            self.p.kill()
            self.join()
            return -9
        return self.retval

class JavaProgram(object):

    def __init__(self,verbose=False,wait=True,stdout=False,javaw=(sys.platform=="win32"),timeout=None):
        self.verbose = verbose
        self.wait = wait
	self.java = 'javaw' if javaw else 'java'
	self.stdout = stdout
        self.timeout = timeout
    
    def __call__(self):
	if 'JAVA_HOME' in os.environ:
	    prefix = os.environ['JAVA_HOME']
	    java = os.path.join(prefix,'bin',self.java)
	else:
	    java = self.java
        cmd = '"%s" -cp "%s" %s %s'%(java,self.classpath(),self.main,self.args())
        if self.verbose:
            print >>sys.stderr, "Executing:", cmd
        starttime = time.time()
        if self.stdout:
            proc = PopenTimeout(cmd,
                                stdin=PIPE,stdout=PIPE,stderr=STDOUT,
                                shell=(sys.platform!="win32"),
                                timeout=self.timeout)
        else:
            proc = PopenTimeout(cmd,
                                stdin=PIPE,
                                shell=(sys.platform!="win32"),
                                timeout=self.timeout)
        proc.p.stdin.write(self.stdin())
        if self.wait == False:
            return
        if type(self.wait) in (int,float):
            time.sleep(self.wait)
            return
        retval = proc.wait()
        if retval == -9:
            if self.verbose:
                print >>sys.stderr, "Process killed after %s seconds"%(self.timeout,)
        else:
            if self.verbose:
                print >>sys.stderr, "Process completed after %s seconds"%(round(time.time()-starttime,1),)
	if self.stdout:
            return proc.p.stdout.read()
	return (retval==0)

    def stdin(self):
        return ""

    def classpath(self):
        classpath = self.findlibs()
        cpsep = ':' if (sys.platform!="win32") else ';'
        classpath = cpsep.join(classpath)
        return classpath

    def findlibs(self):
        bindir = __file__

        while not os.path.isdir(bindir):
            bindir,dummy = os.path.split(bindir)

        found_libs = []
        for l in map(lambda l: l+".jar",self.libs.split()):
            curdir = bindir
            l1 = os.path.join(curdir,l)
            if os.path.exists(l1):
                found_libs.append(l1)
                continue
            curdir = os.path.join(curdir,'..')
            l1 = os.path.join(curdir,l)
            if os.path.exists(l1):
                found_libs.append(l1)
                continue
            curdir = os.path.join(curdir,'lib')
            l1 = os.path.join(curdir,l)
            if os.path.exists(l1):
                found_libs.append(l1)
                continue
            raise RuntimeError("Can't find library %s"%l)

        return found_libs
        
class GlycoCT2Image(JavaProgram):
    libs = """
              GlycoCT2ImageBundle
           """
    main = "GlycoCT2Image"

    def __init__(self,glycoctstr,outfile,verbose=False,timeout=15,**kw):
        super(GlycoCT2Image,self).__init__(verbose=verbose,wait=True,stdout=(not verbose),timeout=timeout)
        self.kwargs = kw
        self.glycoctstr = glycoctstr
        self.outfile = outfile

    def args(self):
        theargs = []
        for k,v in self.kwargs.items():
            theargs.append(k)
            theargs.append(v)
        theargs.append("out")
        theargs.append(self.outfile)
        theargs.append("-")
        return " ".join(map(str,theargs))

    def stdin(self):
        return self.glycoctstr

class GWBFormatter(JavaProgram):
    libs = """
              GlycoCT2ImageBundle
           """
    main = "GlycoCT2GWB"

    def __init__(self,glycoctstr,verbose=False,**kw):
        super(GWBFormatter,self).__init__(verbose=verbose,wait=True,stdout=True,timeout=10)
        self.glycoctstr = glycoctstr

    def args(self):
        return "-"

    def stdin(self):
        return self.glycoctstr

