from __future__ import print_function

from .WebServiceResource import WebServiceResource

import os, os.path, re, time, sys
import hashlib 

class GlycamWS(WebServiceResource):
    apiurl = 'https://glycam.org/json'

    def __init__(self,**kw):
        kw['iniFile'] = os.path.join(os.path.dirname(os.path.realpath(__file__)),"glycam.ini")
        kw['delaytime'] = 0
        super(GlycamWS,self).__init__(**kw)
        token = self.opener.open(self.apiurl+'/getToken/').read().decode()
        self.opener.addheaders.append(('X-CSRFToken', token))
        if self._verbose:
            for c in self.cj:
                print("Cookie %s: %s (%s)"%(c.name,c.value,c.domain),file=sys.stderr)
            print("Added X-CSRFToken: %s"%(token,),file=sys.stderr) 

    def valid(self,*seq):
        for s in seq:
            data = self.query_valid(sequence=s)
            try:
                if data["entity"]["outputs"] == None:
                    yield dict(valid=False)
                else:
                    yield dict(valid=data["entity"]["outputs"]["sequenceEvaluationOutput"]["sequenceIsValid"])
            except (TypeError,KeyError):
                yield dict(valid=None)
