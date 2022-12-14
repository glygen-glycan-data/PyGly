

from .WebServiceResource import WebServiceResource

import os, os.path, re, time
import hashlib 

class GlycoTreeSandbox(WebServiceResource):
    apiurl = 'https://glygen.ccrc.uga.edu/sandbox/api/'

    def __init__(self,**kw):
        kw['iniFile'] = os.path.join(os.path.dirname(os.path.realpath(__file__)),"glycotreews.ini")
        super(GlycoTreeSandbox,self).__init__(**kw)

    def list(self,mode='all'):
        assert mode in ('all','all_N','all_O','mapped_N','mapped_O')
        for d in self.query_list(mode=mode):
            yield d['glytoucan_ac']
