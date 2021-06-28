
from tempfile import mkstemp
import os,sys,re,os.path,time
from subprocess import Popen, PIPE, STDOUT

class JavaProgram(object):

    def __init__(self,verbose=False,wait=True,stdout=False,javaw=(sys.platform=="win32")):
        self.verbose = verbose
        self.wait = wait
	self.java = 'javaw' if javaw else 'java'
	self.stdout = stdout
    
    def __call__(self):
	if 'JAVA_HOME' in os.environ:
	    prefix = os.environ['JAVA_HOME']
	    java = os.path.join(prefix,'bin',self.java)
	else:
	    java = self.java
        cmd = '"%s" -cp "%s" %s %s'%(java,self.classpath(),self.main,self.args())
        if self.verbose:
            print >>sys.stderr, "Executing:", cmd
        if self.stdout:
            proc = Popen(cmd,
                         stdin=PIPE,stdout=PIPE,stderr=STDOUT,
                         shell=(sys.platform!="win32"))
        else:
            proc = Popen(cmd,
                         stdin=PIPE,
                         shell=(sys.platform!="win32"))
        proc.stdin.write(self.stdin())
        if self.wait == False:
            return
        if type(self.wait) in (int,float):
            time.sleep(self.wait)
            return
        retval = proc.wait()
	if self.stdout:
            return proc.stdout.read()
	return (retval!=0)

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

    def __init__(self,glycoctstr,outfile,verbose=False,wait=True,stdout=False,**kw):
        super(GlycoCT2Image,self).__init__(verbose=verbose,wait=wait,stdout=stdout)
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

    def __init__(self,glycoctstr,verbose=False,wait=True,**kw):
        super(GWBFormatter,self).__init__(verbose=verbose,wait=wait,stdout=True)
        self.glycoctstr = glycoctstr

    def args(self):
        return "-"

    def stdin(self):
        return self.glycoctstr

