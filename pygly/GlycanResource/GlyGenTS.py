
import os

from TripleStoreResource import TripleStoreResource
from GlycanResourceWrappers import partitioner

class GlyGenTS(TripleStoreResource):
    endpt = "http://sparql.glygen.org:8880/sparql/query"
    defns = "http://glygen.org/glycan/"

    def __init__(self,**kw):
        kw['iniFile'] = os.path.join(os.path.dirname(os.path.realpath(__file__)),"glygen.ini")
        super(GlyGenTS,self).__init__(**kw)
        for k in self.keys():
            self.modify_method(k,partitioner())

    def allglycans(self):
	for row in self.query_glycans():
            yield row['accession']

class GlyGenBetaTS(GlyGenTS):
    endpt = "http://beta-sparql.glygen.org:8880/sparql/query"
    def __init__(self):
        super(GlyGenBetaTS,self).__init__(usecache=False)
