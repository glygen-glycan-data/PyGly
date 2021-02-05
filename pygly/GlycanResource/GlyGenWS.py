
from WebServiceResource import WebServiceResource

import os, os.path, re, time
import hashlib 

class GlyGenWS(WebServiceResource):
    apiurl = 'https://api.glygen.org/'
		
    def __init__(self,**kw):
	kw['iniFile'] = os.path.join(os.path.dirname(os.path.realpath(__file__)),"glygenws.ini")
	super(GlyGenWS,self).__init__(**kw)
