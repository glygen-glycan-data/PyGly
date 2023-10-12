
from .GlycoMotifTS import GlycoMotifTS, GlycoMotifDevTS

class GlycoMotif(GlycoMotifTS):
    pass

class GlycoMotifNoPrefetch(GlycoMotif):
    def __init__(self,**kw):
        kw['prefetch'] = False
        super(GlycoMotifNoPrefetch,self).__init__(**kw);

class GlycoMotifNoCache(GlycoMotif):
    def __init__(self,**kw):
        kw['usecache'] = False
        super(GlycoMotifNoCache,self).__init__(**kw);

class GlycoMotifDev(GlycoMotifDevTS):
    pass

class GlycoMotifDevNoPrefetch(GlycoMotifDev):
    def __init__(self,**kw):
        kw['prefetch'] = False
        super(GlycoMotifDevNoPrefetch,self).__init__(**kw);

class GlycoMotifDevNoCache(GlycoMotifDev):
    def __init__(self,**kw):
        kw['usecache'] = False
        super(GlycoMotifDevNoCache,self).__init__(**kw);

