
from __future__ import print_function
from .GlycanResource import GlycanResource


import os, ssl

if (not os.environ.get('PYTHONHTTPSVERIFY', '') and getattr(ssl, '_create_unverified_context', None)):
    ssl._create_default_https_context = ssl._create_unverified_context

import sys, json, csv

try:
    from urllib.parse import urlparse, urlencode
    from urllib.request import urlopen, Request, build_opener, HTTPSHandler, HTTPHandler, HTTPCookieProcessor
    from urllib.error import HTTPError
    from http.cookiejar import CookieJar
except ImportError:
    from urlparse import urlparse
    from urllib import urlencode
    from urllib2 import urlopen, Request, HTTPError, build_opener, HTTPSHandler, HTTPHandler, HTTPCookieProcessor
    from cookielib import CookieJar

from io import StringIO, TextIOWrapper, BytesIO

try:
    from gzip import decompress as gzipdecompress
except ImportError:
    import zlib
    def gzipdecompress(b):
        return zlib.decompress(b,16+zlib.MAX_WBITS)

class WebServiceResource(GlycanResource):

    """
    Class for accessing glycan data from webservices

    url: Base URL for webservices

    """

    def __init__(self,*args,**kw):
        self.attr(kw,'apiurl',required=True)
        self.attr(kw,'verbose',default=False)
        self.attr(kw,'local',default=False)
        self.attr(kw,'localurl',default=None)

        if self._local:
            self._apiurl = self._localurl

        self.cj = CookieJar()
        self.opener = build_opener(HTTPCookieProcessor(self.cj))
        # self.opener = build_opener(HTTPCookieProcessor(self.cj),HTTPSHandler(debuglevel=1))
        super(WebServiceResource,self).__init__(*args,**kw)

    def queryws(self,url,method,kwargs):
        self.wait()
        if method == "GET":
          if len(kwargs) > 0:
            if '?' in url or '%(' in url:        
                url = url%kwargs
            else:
                url += "?" + urlencode(kwargs)
          if self._verbose:
            print(url,file=sys.stderr)
          req = Request(url)
          h = self.opener.open(req)
        elif method == "POST":
          payload = kwargs.get('payload','{}')
          payload = payload%kwargs
          payload = json.loads(payload)
          if self._verbose:
            print(url,file=sys.stderr)
            print(json.dumps(payload),file=sys.stderr)
          req = Request(url,data=json.dumps(payload).encode(),headers={'Content-Type': 'application/json'})
          h = self.opener.open(req)
        return h.read()

    def parseSection(self,name,keyvaluepairs):
        method = keyvaluepairs.get("method","GET")
        url = keyvaluepairs.get("url")
        if not url:
            req = keyvaluepairs.get("request",name)
            url = self._apiurl + "/" + req
        params = [_f for _f in [s.strip() for s in keyvaluepairs.get('params',"").split(',')] if _f] + ["payload"]
        thetype = keyvaluepairs.get("type","JSON").upper()

        def _query(self,*args,**kw):
            kwargs = {}
            for i,param in enumerate(params):
                if param in keyvaluepairs:
                    kwargs[param] = keyvaluepairs[param]
                if i < len(args):
                    kwargs[param] = args[i]
                elif kw.get(param) != None:                                                                                 
                    kwargs[param] = kw[param]                                                                               
                assert param in kwargs or param == "payload", " ".join(map(repr,[param, kwargs])) 
            response = self.queryws(url,method,kwargs)
            if thetype == "JSON":
                response = json.loads(response)
            elif thetype == "CSV":
                response = csv.DictReader(StringIO(response.decode(encoding='ascii',errors='ignore')))
            elif thetype == "TSV":
                response = csv.DictReader(StringIO(response.decode(encoding='ascii',errors='ignore')),dialect='excel-tab')
            elif thetype == "TEXTGZ":
                response = TextIOWrapper(BytesIO(gzipdecompress(response)),encoding='utf8')
            elif thetype == "TEXT":
                try:
                    response = TextIOWrapper(response,encoding='utf8')
                except AttributeError:
                    pass
            return response

        self.set_method("query_"+str(name), _query)
        return [("query_"+name,params)]
