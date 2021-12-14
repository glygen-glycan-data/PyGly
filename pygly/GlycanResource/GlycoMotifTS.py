
import os.path
from .TripleStoreResource import TripleStoreResource
from .GlycanResourceWrappers import partitioner, prefetcher

class GlycoMotifTS(TripleStoreResource):
   
    endpt = "https://glycomotif.glyomics.org/glycomotif/sparql/query"
    defns = "http://glyomics.org/glycomotif/Special:URIResolver/"
    # verbose = True
    # cachefile = ".gm.cache"
    # usecache = False
    # cachemode= 'c'
    def __init__(self,**kw):
        if 'prefetch' not in kw:
            kw['prefetch'] = True
        kw['iniFile'] = os.path.join(os.path.dirname(os.path.realpath(__file__)),"glycomotif.ini")
        super(GlycoMotifTS,self).__init__(**kw)
        for k in list(self.keys()):
            if k != 'query_motifs':
                continue
            self.modify_method(k,partitioner())
            if self._prefetch:
                self.modify_method(k,prefetcher(usecache=False))

    def getmotif(self,collection,accession):
        return [ ("%s.%s"%(collection,row['MotifAccession']),row['StrictAlignment']=='true') for row in self.query_motifs(collection=collection,accession=accession) ]

    def allmotifaligns(self,collection):
        for row in self.query_motifs(collection=collection):
            yield row['accession'],"%s.%s"%(collection,row['MotifAccession']),row['StrictAlignment']=='true'

    def allmotifs(self,collection):
        for row in self.query_allmotif(collection=collection):
            names = []
            if row.get('prefname'):
                names.append(row.get('prefname'))
            if row.get('name'):
                for name in row.get('name').split('//'):
                    if name not in names:
                        names.append(name)
            pmids = []
            if row.get('pmid'):
                pmids = row.get('pmid').split('//')
            keywords = []
            if row.get('keyword'):
                keywords = row.get('keyword').split('//')
            dbxrefs = []
            if row.get('dbxref'):
                dbxrefs = row.get('dbxref').split('//')
            yield "%s.%s"%(collection,row['accession']),row['gtcacc'],row['alignment'],row['redend'],row['aglycon'],names,pmids,keywords,dbxrefs
