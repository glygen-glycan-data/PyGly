
import os.path
from .TripleStoreResource import TripleStoreResource
from .GlycanResourceWrappers import partitioner, prefetcher

class GlycoMotifTS(TripleStoreResource):
   
    endpt = "https://glycomotif.glyomics.org/glycomotif/sparql/query"
    localendpt = "http://trypsin:8081/glycomotif/sparql/query"
    defns = "http://glyomics.org/glycomotif/Special:URIResolver/"
    # verbose = True
    cachefile = ".gm.cache"
    # usecache = True
    # cachemode= 'c'
    def __init__(self,**kw):
        if 'prefetch' not in kw:
            kw['prefetch'] = True
        if 'usecache' not in kw:
            kw['usecache'] = True
        if 'cachemode' not in kw:
            kw['cachemode'] = 'c'
        kw['iniFile'] = os.path.join(os.path.dirname(os.path.realpath(__file__)),"glycomotif.ini")
        kw['delaytime'] = 0.5
        kw['delaybatch'] = 2
        super(GlycoMotifTS,self).__init__(**kw)
        for k in list(self.keys()):
            if k == 'query_motifs':
                self.modify_method(k,partitioner())
            elif k == 'query_struct':
                self.modify_method(k,partitioner(kwarg="motifacc",fmt=".*%%0%dd"))
            elif k == 'query_hasstruct':
                self.modify_method(k,partitioner())
            if self._prefetch:
                if k == 'query_motifs':
                    self.modify_method(k,prefetcher(usecache=self._usecache))
                elif k == 'query_struct':
                    self.modify_method(k,prefetcher("motifacc",usecache=self._usecache))
                elif k == 'query_hasstruct':
                    self.modify_method(k,prefetcher(usecache=self._usecache))

    def collections(self):
        for row in self.query_collections():
            yield row['CollectionID']

    def hasstructure(self,accession):
        for row in self.query_hasstruct(accession=accession):
            return True
        return False

    def getmotif(self,collection,accession):
        return [ ("%s.%s"%(collection,row['MotifAccession']),row['MotifAlignment'],row['StrictAlignment']=='true',row['StructureResidueIds'],row['StructureLinkIds']) for row in self.query_motifs(collection=collection,accession=accession) ]

    def getstruct(self,collection,motifacc):
        return sorted(set((row['accession'],row['StrictAlignment']=='true',row['StructureResidueIds'],row['StructureLinkIds']) for row in self.query_struct(collection=collection,motifacc=motifacc)))

    def getmotifbystruct(self,accession):
        return [ ("%s.%s"%(row['CollectionID'],row['MotifAccession']),row['MotifAlignment'],row['StrictAlignment']=='true',row['StructureResidueIds'],row['StructureLinkIds']) for row in self.query_motifsbystructure(accession=accession) ]

    def allmotifaligns(self,collection):
        for row in self.query_motifs(collection=collection):
            yield row['accession'],"%s.%s"%(collection,row['MotifAccession']),row['MotifAlignment'],row['StrictAlignment']=='true',row['StructureResidueIds'],row['StructureLinkIds']

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

class GlycoMotifDevTS(GlycoMotifTS):
   
    endpt = "https://glycomotif.glyomics.org/glycomotifdev/sparql/query"
    localendpt = "http://trypsin:8080/glycomotifdev/sparql/query"
    defns = "http://glyomics.org/glycomotifdev/Special:URIResolver/"
    cachefile = ".gmd.cache"

