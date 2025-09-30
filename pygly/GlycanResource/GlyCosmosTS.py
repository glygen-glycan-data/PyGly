
from .TripleStoreResource import TripleStoreResource
from .GlycanResourceWrappers import partitioner, prefetcher

import hashlib
import base64
from collections import defaultdict

class GlyCosmosTS(TripleStoreResource):

    endpt = "https://ts.glycosmos.org/sparql"
    # endpt = "http://ts.beta.glycosmos.org/sparql"
    defns = "http://rdf.glycoinfo.org/glycan/"
    # verbose = True
    cachefile = ".gco.cache"
    # usecache = True
    # prefetch = True
    # cachemode= 'c'

    sequence_formats = ['wurcs','glycoct','iupac_extended','iupac_condensed']
    
    def __init__(self,**kw):
        if 'usecache' not in kw:
            kw['usecache'] = True
        if 'prefetch' not in kw:
            kw['prefetch'] = True
        if 'cachemode' not in kw:
            kw['cachemode'] = 'c'
        super(GlyCosmosTS,self).__init__(**kw)
        for k in list(self.keys()):
            if 'hash' in self[k]:
                self.modify_method(k,partitioner(kwarg="hash",fmt="%%0%dx.*",values='hexidecimal'))
            elif 'accession' in self[k]:
                self.modify_method(k,partitioner())
            if self._prefetch:
                if 'hash' in self[k]:
                    self.modify_method(k,prefetcher("hash",usecache=self._usecache))
                elif 'accession' in self[k]:
                    self.modify_method(k,prefetcher(usecache=self._usecache))

    def exists(self,accession):
        for row in self.query_exists(accession=accession):
            return True
        return False

    def allaccessions(self):
        for row in self.query_exists():
            yield row['accession']

    def getseq(self,accession,format='wurcs'):
        assert format in self.sequence_formats
        for row in self.query_sequence(accession=accession,format=format):
            if row['format'] == 'wurcs' and not row['sequence'].startswith('WURCS'):
                continue
            if row['format'] == 'glycoct' and not row['sequence'].startswith('RES'):
                continue
            return row['sequence']
        return None

    def allseq(self,format=None):
        assert format == None or format in self.sequence_formats
        curkey = None
        for row in self.query_sequence(format=format):
            if row['format'] == 'wurcs' and not row['sequence'].startswith('WURCS'):
                continue
            if row['format'] == 'glycoct' and not row['sequence'].startswith('RES'):
                continue
            key = (row['accession'],row['format'])
            if key != curkey:
                yield row['accession'], row['format'], row['sequence']
                curkey = key

    def getmass(self,accession):
        for row in self.query_mass(accession=accession):
            try:
                return float(row['mass'])
            except ValueError:
                pass
        return None

    def allmass(self):
        for row in self.query_mass():
            try:
                yield row['accession'],float(row['mass'])
            except ValueError:
                pass

    taxa_sources = {
       'DOI': True,
       'PubMed': True,

       'GlycomeAtlas': True,

       'GlyTouCan(glycome-db)': True,
       'GlyTouCan(bcsdb)': True,

       'GPTwiki': True,

       'GlycoEpitope': True,

       'GlyConnectStructure': False,
       'GlyConnectComposition': False,

       'GlyGen(UniCarbKB)': False,
       'GlyGen(OGlcNAcAtlas)': False,
       'GlyGen(OGlcNAcDB)': False,
       'GlyGen(Harvard)': False,
       'GlyGen': False,
    }

    def gettaxa(self,accession):
        for row in self.query_taxonomy(accession=accession):
            if self.taxa_sources[row['source']]:
                yield row['taxon']

    def bytaxa(self,taxon,*extrataxon):
        if len(extrataxon) > 0:
            taxonre = '(' + "|".join(map(str,[ taxon ] + list(extrataxon))) + ')'
        else:
            taxonre = str(taxon)
        for row in self.query_taxonomy(taxon=taxonre):
            if self.taxa_sources[row['source']]:
                yield row['accession'],row['taxon']

    def alltaxa(self):
        for row in self.query_taxonomy():
            if self.taxa_sources[row['source']]:
                yield row['accession'],row['taxon']

