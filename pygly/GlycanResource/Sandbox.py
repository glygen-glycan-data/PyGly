

from .WebServiceResource import WebServiceResource

import os, os.path, re, time, json
import hashlib 

class GlycoTreeSandbox(WebServiceResource):
    # apiurl = 'https://glygen.ccrc.uga.edu/sandbox/api/'
    apiurl = 'https://sandbox.glyomics.org/api'

    def __init__(self,**kw):
        kw['iniFile'] = os.path.join(os.path.dirname(os.path.realpath(__file__)),"glycotreews.ini")
        super(GlycoTreeSandbox,self).__init__(**kw)

    def list(self,mode='all'):
        assert mode in ('all','all_N','all_O','mapped','mapped_N','mapped_O','clean','clean_N','clean_O')
        for d in self.query_list(mode=mode):
            yield d['glytoucan_ac']

    def glycan(self,accession):
        g = self.query_glycan(accession=accession)
        g['accession'] = g['glytoucan_ac']
        del g['glytoucan_ac']
        return g

    def glycans(self,*accs):
        for g in self.query_glycans(accessions=",".join(accs)):
            g['accession'] = g['glytoucan_ac']
            del g['glytoucan_ac']
            yield g

    def allglycans(self,mode='all',blocksize=20):
        assert mode in ('all','all_N','all_O','mapped','mapped_N','mapped_O','clean','clean_N','clean_O')
        listaccs = list(self.list(mode))
        for i in range(0,len(listaccs),blocksize):
            accs = listaccs[i:(i+blocksize)]
            for r in self.glycans(*accs):
                yield r

class GlycoTreeSandboxDev(GlycoTreeSandbox):
    apiurl = 'https://edwardslab.bmcb.georgetown.edu/sandboxdev/api'

