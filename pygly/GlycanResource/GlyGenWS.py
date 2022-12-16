
from .WebServiceResource import WebServiceResource

import os, os.path, re, time, json
import hashlib 

class GlyGenWS(WebServiceResource):
    apiurl = 'https://api.glygen.org/'
                
    def __init__(self,**kw):
        kw['iniFile'] = os.path.join(os.path.dirname(os.path.realpath(__file__)),"glygenws.ini")
        super(GlyGenWS,self).__init__(**kw)

    def glycan_search(self,**kw):
        offset = 1
        while True:
            query=json.dumps(dict(offset=offset,**kw))
            offset += 1
            result = self.query_glycan_search(query=query)
            print(dict(len=len(result['results']),**result['queryinfo']))
            for d in result['results']:
                yield d
            if len(result['results']) == 0:
                break

    def glycans_bytype(self,type):
         for i,d in enumerate(self.glycan_search(glycan_type=type)):
             yield d['glytoucan_ac']
