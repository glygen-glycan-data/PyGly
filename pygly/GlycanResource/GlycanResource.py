
import time
try:
    from pygly.ReferenceTable import ReferenceTable
except ImportError:
    from ReferenceTable import ReferenceTable

class GlycanResource(ReferenceTable):
    """
       Abstract base class for glycan resources whose data is
       accessed using triple store queries or HTTP requests

       delaytime: Waiting time between request batches, default 0.2 sec
       delaybatch: Size of request batches with no waiting time, default 1 (no batches)
       retries: Maximum number of retries: default 4

    """

    def __init__(self,**kw):
        self._delaybatch = kw.get('delaybatch',1)
        self._delaytime = kw.get('delaytime',0.2)
        self._retries = kw.get('retries',4)
        self._lastrequesttime = 0
        self._requestcount = 0

        super(GlycanResource,self).__init__(iniFile=kw.get('iniFile'))

    def wait(self,delay=None):
        elapsed = time.time() - self._lastrequesttime
        if delay != None:
            if elapsed < delay:
                time.sleep(delay-elapsed)
        elif (self._requestcount % self._delaybatch) == 0 and self._requestcount > 0:
            if elapsed < delay:
                time.sleep(self._delaytime-elapsed)
        self._lastrequesttime = time.time()
        self._requestcount += 1

    def attr(self,kw,key,default=None,required=False):
        if hasattr(self,key):
            setattr(self,"_"+key,getattr(self,key))
        elif key in kw:
            setattr(self,"_"+key,kw[key])
        elif not required:
            setattr(self,"_"+key,default)
        else:
            raise RuntimeError("Can't find class/instance parameter %s for class %s"%(key,self.__class__.__name__))

    def set_method(self,name,func):
        setattr(self.__class__, name, func)
        func.__name__ = name
    
    def modify_method(self,name,func):
        newfunc = func(getattr(self.__class__,name))
        setattr(self.__class__, name, newfunc)
        newfunc.__name__ = name
