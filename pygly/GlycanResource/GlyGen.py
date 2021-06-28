
from .GlyGenTS import GlyGenTS, GlyGenBetaTS
from .GlyGenWS import GlyGenWS

class GlyGen(GlyGenTS):
    def __init__(self,**kw):
        super(GlyGen,self).__init__(**kw)
        self.ws = GlyGenWS()

    def alltimeglycans(self):
        for acc in self.ws.query_alltime():
            yield acc

class GlyGenBeta(GlyGenBetaTS):
    pass
