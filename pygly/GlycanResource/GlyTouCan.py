
from .GlyTouCanTS import GlyTouCanTS
from .GlyTouCanUtil import GlyTouCanUtil
from .GlyTouCanRegistration import GlyTouCanRegistration
from .GlycanResourceWrappers import cacher

import os.path

class GlyTouCan(GlyTouCanTS,GlyTouCanUtil,GlyTouCanRegistration):
    def __init__(self,**kw):
        kw['iniFile'] = os.path.join(os.path.dirname(os.path.realpath(__file__)),"glytoucan.ini")
        super(GlyTouCan,self).__init__(**kw)
        self.modify_method("umw",cacher(usecache=self._usecache))

class GlyTouCanNoCache(GlyTouCan):
    def __init__(self,**kw):
        kw['usecache'] = False
        super(GlyTouCanNoCache,self).__init__(**kw);

class GlyTouCanNoPrefetch(GlyTouCan):
    def __init__(self,**kw):
        kw['prefetch'] = False
        super(GlyTouCanNoPrefetch,self).__init__(**kw);
