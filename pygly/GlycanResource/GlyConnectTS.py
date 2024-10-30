
from .TripleStoreResource import TripleStoreResource
from .GlycanResourceWrappers import partitioner, prefetcher

import hashlib, os, os.path
import base64
from collections import defaultdict

class GlyConnectTS(TripleStoreResource):

    endpt = "https://glyconnect.expasy.org/sparql"
    defns = "http://identifiers.org/glytoucan/"
    method = "POST"
    # verbose = True
    cachefile = ".gcn.cache"
    # usecache = True
    # prefetch = True
    # cachemode= 'c'

    def __init__(self,**kw):
        if 'usecache' not in kw:
            kw['usecache'] = True
        if 'prefetch' not in kw:
            kw['prefetch'] = True
        if 'cachemode' not in kw:
            kw['cachemode'] = 'c'
        kw['iniFile'] = os.path.join(os.path.dirname(os.path.realpath(__file__)),"glyconnectts.ini")
        super(GlyConnectTS,self).__init__(**kw)
        for k in list(self.keys()):
            if 'accession' in self[k]:
                self.modify_method(k,partitioner())
            if self._prefetch:
                if 'accession' in self[k]:
                    self.modify_method(k,prefetcher(usecache=self._usecache))

    def allaccessions(self):
        for row in self.query_exists():
            yield row['accession']

    def exists(self,accession):
        for row in GlyConnectTS.query_exists(self,accession=accession):
            return True
        return False

    def allgtc(self):
        for row in self.query_exists():
            yield row['gcnid'],row['accession']

class GlyConnectTSNoCache(GlyConnectTS):
    def __init__(self,**kw):
        kw['usecache'] = False
        GlyConnectTS.__init__(self,**kw)

class GlyConnectTSNoPrefetch(GlyConnectTS):
    def __init__(self,**kw):
        kw['prefetch'] = False
        GlyConnectTS.__init__(self,**kw)
