

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
        assert mode in ('all','all_N','all_O','mapped_N','mapped_O')
        for d in self.query_list(mode=mode):
            yield d['glytoucan_ac']

    def glycan(self,accession):
        g = self.query_glycan(accession=accession)
        g['accession'] = g['glytoucan_ac']
        del g['glytoucan_ac']
        return g

    def allglycans(self,mode='all'):
        assert mode in ('all','all_N','all_O','mapped_N','mapped_O')
        for acc in self.list(mode):
            yield self.glycan(acc)

class GlycoTreeSandboxDev(GlycoTreeSandbox):
    apiurl = 'https://edwardslab.bmcb.georgetown.edu/sandboxdev/api'
