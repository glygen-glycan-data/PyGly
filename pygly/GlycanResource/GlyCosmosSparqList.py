
from .WebServiceResource import WebServiceResource

import os, os.path, re, time
import hashlib 

class GlyCosmosSparqList(WebServiceResource):
    apiurl = 'https://sparqlist.glycosmos.org/sparqlist/api'
    export = ['archived','replaced','validated','replace']

    def __init__(self,**kw):
        kw['iniFile'] = os.path.join(os.path.dirname(os.path.realpath(__file__)),"glycosmos_sparqlist.ini")
        super(GlyCosmosSparqList,self).__init__(**kw)

    def archived(self):
        for it in self.query_archived():
            yield dict(accession=it['ArchiveNumber'])

    def replaced(self):
        for it in self.query_replaced():
            yield dict(accession=it['AccessionNumber'],replace=it['ArchiveNumber'])

    def validated(self):
        seen = set()
        for it in self.query_validated():
            if it['glycan'] in seen:
                continue
            seen.add(it['glycan'])
            yield dict(accession=it['glycan'])

    def replace(self):
        repl = dict()
        for it in self.query_archived():
            repl[it['ArchiveNumber']] = None
        for it in self.query_replaced():
            repl[it['ArchiveNumber']] = it['AccessionNumber']
        return repl

