
from .GlyCosmosTS import GlyCosmosTS
from .GlyTouCanUtil import GlyTouCanUtil
from .GlyCosmosSparqList import GlyCosmosSparqList

import os.path

class GlyCosmos(GlyCosmosTS,GlyTouCanUtil):
    def __init__(self,**kw):
        kw['iniFile'] = os.path.join(os.path.dirname(os.path.realpath(__file__)),"glycosmos.ini")
        super(GlyCosmos,self).__init__(**kw)
        self.sparqlist = GlyCosmosSparqList()
        self.attach_methods(self.sparqlist)

class GlyCosmosNoCache(GlyCosmos):
    def __init__(self,**kw):
        kw['usecache'] = False
        super(GlyCosmosNoCache,self).__init__(**kw);

class GlyCosmosNoPrefetch(GlyCosmos):
    def __init__(self,**kw):
        kw['prefetch'] = False
        super(GlyCosmosNoPrefetch,self).__init__(**kw);
