
from .WebServiceResource import WebServiceResource

import os, os.path, re, time, json
import hashlib 

class GlyGenWS(WebServiceResource):
    apiurl = 'https://api.glygen.org'
                
    def __init__(self,**kw):
        kw['iniFile'] = os.path.join(os.path.dirname(os.path.realpath(__file__)),"glygenws.ini")
        super(GlyGenWS,self).__init__(**kw)

    def dataset(self, filename, version=None):
        for row in self.query_csvdataset(filename=filename, version=version):
            yield row

    def protein_homolog_clusters(self, version=None):
        for row in self.dataset("protein_homolog_clusters",version=version):
            yield row

    def protein_genenames(self, species, version=None):
        for row in self.dataset(species+"_protein_genenames_uniprotkb",version=version):
            yield row

    def glycosyltransferases(self, species, version=None):
        for row in self.dataset(species+"_protein_glycosyltransferase",version=version):
            yield row

    def protein_geneid(self, species, version=None):
        for row in self.dataset(species+"_protein_xref_geneid",version=version):
            yield row

    def protein_refseqnp(self, species, version=None):
        for row in self.dataset(species+"_protein_xref_refseq",version=version):
            yield row

    def protein_masterlist(self, species, version=None):
        for row in self.dataset(species+"_protein_masterlist",version=version):
            yield row

    def glycan_search(self, **kw):
        query=json.dumps(kw)
        result = self.query_glycan_search(query=query)
        list_id = result['list_id']
        result = self.query_get_list(search_id=list_id)
        listcache_id=result['cache_info']['listcache_id']
        headers = None
        for l in self.query_download_list(list_id=listcache_id):
            if not headers:
                headers = list(map(lambda s: s.strip('"'),l.split(",")))
                continue
            sl = l.strip('"').split('","')
            yield dict(zip(headers,sl))

    def glycan_directsearch(self,**kw):
        offset = 1
        while True:
            query=json.dumps(dict(offset=offset,**kw))
            offset += 1
            result = self.query_glycan_directsearch(query=query)
            for d in result['results']:
                yield d
            if len(result['results']) == 0:
                break

    def glycans_bytype(self,type):
         for i,d in enumerate(self.glycan_search(glycan_type=type)):
             yield d['glytoucan_ac']
