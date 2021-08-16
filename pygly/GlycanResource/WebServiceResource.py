
from .GlycanResource import GlycanResource


import os, ssl

if (not os.environ.get('PYTHONHTTPSVERIFY', '') and getattr(ssl, '_create_unverified_context', None)):
    ssl._create_default_https_context = ssl._create_unverified_context

import sys, json, csv

try:
    from urllib.parse import urlparse, urlencode
    from urllib.request import urlopen, Request, build_opener, HTTPSHandler, HTTPHandler
    from urllib.error import HTTPError
except ImportError:
    from urlparse import urlparse
    from urllib import urlencode
    from urllib2 import urlopen, Request, HTTPError, build_opener, HTTPSHandler, HTTPHandler

from io import StringIO

class WebServiceResource(GlycanResource):

    """
    Class for accessing glycan data from webservices

    url: Base URL for webservices

    """

    def __init__(self,*args,**kw):
        self.attr(kw,'apiurl',required=True)
        self.attr(kw,'verbose',default=False)
        super(WebServiceResource,self).__init__(*args,**kw)

    def queryws(self,url,method,kwargs):
        self.wait()
        assert method == "GET"
        if len(kwargs) > 0:
            if '?' in url or '%(' in url:        
                url = url%kwargs
            else:
                url += "?" + urlencode(kwargs)
        if self._verbose:
            print >>sys.stderr, url
        # req = Request(url)
        # req.add_header('User-Agent', 'Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:88.0) Gecko/20100101 Firefox/88.0')
        h = urlopen(url)
        return h.read()

    def parseSection(self,name,keyvaluepairs):
        method = keyvaluepairs.get("method","GET")
        url = keyvaluepairs.get("url")
        if not url:
            req = keyvaluepairs.get("request",name)
            url = self._apiurl + "/" + req
        params = [_f for _f in [s.strip() for s in keyvaluepairs.get('params',"").split(',')] if _f]
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
                assert param in kwargs, " ".join(map(repr,[param, kwargs]))                                                 
            response = self.queryws(url,method,kwargs)
            if thetype == "JSON":
                response = json.loads(response)
            elif thetype == "CSV":
                response = csv.DictReader(StringIO(response.decode(encoding='ascii',errors='ignore')))
            return response

        self.set_method("query_"+str(name), _query)
        return [("query_"+name,params)]
