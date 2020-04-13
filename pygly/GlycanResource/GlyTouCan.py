
from GlyTouCanTS import GlyTouCanTS
from GlyTouCanUtil import GlyTouCanUtil
from GlyTouCanRegistration import GlyTouCanRegistration

import os.path

class GlyTouCan(GlyTouCanTS,GlyTouCanUtil,GlyTouCanRegistration):
    def __init__(self,**kw):
	kw['iniFile'] = os.path.join(os.path.dirname(os.path.realpath(__file__)),"glytoucan.ini")
	super(GlyTouCan,self).__init__(**kw)

class GlyTouCanNoCache(GlyTouCan):
    def __init__(self,**kw):
	kw['usecache'] = False
	super(GlyTouCanNoCache,self).__init__(**kw);

class GlyTouCanNoPrefetch(GlyTouCan):
    def __init__(self,**kw):
	kw['prefetch'] = False
	super(GlyTouCanNoPrefetch,self).__init__(**kw);
