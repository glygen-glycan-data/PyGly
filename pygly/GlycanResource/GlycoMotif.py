
from .GlycoMotifTS import GlycoMotifTS

class GlycoMotif(GlycoMotifTS):
    pass

class GlycoMotifNoPrefetch(GlycoMotif):
    def __init__(self,**kw):
        kw['prefetch'] = False
        super(GlycoMotifNoPrefetch,self).__init__(**kw);

