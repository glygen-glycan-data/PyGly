
from .GlyGenTS import GlyGenTS, GlyGenBetaTS
from .GlyGenWS import GlyGenWS

class GlyGen(GlyGenTS):
    def __init__(self,**kw):
        super(GlyGen,self).__init__(**kw)
        self.ws = GlyGenWS(**kw)

    def alltimeglycans(self):
        for acc in self.ws.query_alltime():
            yield acc

    def glycans_bytype(self,type):
         for d in self.ws.glycan_search(glycan_type=type):
             yield d['glytoucan_ac']

class GlyGenBeta(GlyGenBetaTS):
    pass
