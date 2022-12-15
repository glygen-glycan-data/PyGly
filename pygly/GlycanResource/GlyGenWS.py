
from .WebServiceResource import WebServiceResource

import os, os.path, re, time, json
import hashlib 

class GlyGenWS(WebServiceResource):
    apiurl = 'https://api.glygen.org/'
                
    def __init__(self,**kw):
        kw['iniFile'] = os.path.join(os.path.dirname(os.path.realpath(__file__)),"glygenws.ini")
        super(GlyGenWS,self).__init__(**kw)

    def glycan_search(self,**kw):
        query=json.dumps(kw)
        result = self.query_glycan_search(query=query)
        # print(result['queryinfo'])
        return result['results']

    def glycans_bytype(self,type):
         for d in self.glycan_search(glycan_type=type):
             yield d['glytoucan_ac']
