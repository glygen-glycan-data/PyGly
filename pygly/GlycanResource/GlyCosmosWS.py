
from .WebServiceResource import WebServiceResource

import os, os.path, re, time
import hashlib 

class GlyCosmosWS(WebServiceResource):
    apiurl = 'https://glycosmos.org'
    export = ['allaccessions','allmass','getmass','getseq']

    def __init__(self,**kw):
        kw['iniFile'] = os.path.join(os.path.dirname(os.path.realpath(__file__)),"glycosmosws.ini")
        super(GlyCosmosWS,self).__init__(**kw)
        self._cache = None
        self.make_cache()

    def make_cache(self):
        if self._cache == None:
            self._cache = dict()
            for r in self.query_glycans():
                acc = r['Accession Number']
                try:
                    mass = float(r['Monoisotopic mass'])
                except:
                    mass = None
                self._cache[acc] = dict(acc=acc,mass=mass,wurcs=r['WURCS'],**r)

    def allrows(self):
        for d in self._cache.values():
            yield d

    def getrow(self,acc):
        return self._cache.get(acc)

    def allaccessions(self):
        for d in self.allrows():
            yield d['acc']

    def allmass(self):
        for d in self.allrows():
            yield d['acc'],d['mass']

    def allseq(self,format='wurcs'):
        assert(format == 'wurcs')
        for d in self.allrows():
            yield d['acc'],d['WURCS']

    def getmass(self,acc):
        d = self.getrow(acc) 
        if d:
            return d['mass']
        return None

    def getseq(self,acc,format='wurcs'):
        assert(format == 'wurcs')
        d = self.getrow(acc) 
        if d:
            return d['WURCS']
        return None
