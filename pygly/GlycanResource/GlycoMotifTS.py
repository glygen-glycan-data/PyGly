
import os.path
from TripleStoreResource import TripleStoreResource
from GlycanResourceWrappers import partitioner, prefetcher

class GlycoMotifTS(TripleStoreResource):
   
    endpt = "https://edwardslab.bmcb.georgetown.edu/sparql/glycomotif/query"
    defns = "http://glycandata.glygen.org/glycomotif/Special:URIResolver/"
    # verbose = True
    # cachefile = ".gm.cache"
    # usecache = False
    # cachemode= 'c'
    def __init__(self,**kw):
        if 'prefetch' not in kw:
            kw['prefetch'] = True
	kw['iniFile'] = os.path.join(os.path.dirname(os.path.realpath(__file__)),"glycomotif.ini")
	super(GlycoMotifTS,self).__init__(**kw)
        for k in self.keys():
	    if k != 'query_motifs':
		continue
            self.modify_method(k,partitioner())
            if self._prefetch:
                self.modify_method(k,prefetcher(usecache=False))

    def getmotif(self,collection,accession):
        return [ "%s.%s"%(row['MotifCollection'],row['MotifAccession']) for row in self.query_motifs(collection=collection,accession=accession) ]

    def allmotifaligns(self,collection):
        for row in self.query_motifs(collection=collection):
            yield row['accession'],"%s.%s"%(row['MotifCollection'],row['MotifAccession'])

    def allmotifs(self,collection):
        for row in self.query_allmotif(collection=collection):
	    names = []
	    if row.get('preferred_name'):
		names.append(row.get('preferred_name'))
	    if row.get('name'):
		names.extend(row.get('name').split('//'))
            yield "%s.%s"%(collection,row['accession']),row['gtcacc'],row['redend'],row['aglycon'],names
