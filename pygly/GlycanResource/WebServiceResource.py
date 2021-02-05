
from GlycanResource import GlycanResource


import os, ssl
if (not os.environ.get('PYTHONHTTPSVERIFY', '') and getattr(ssl, '_create_unverified_context', None)):
    ssl._create_default_https_context = ssl._create_unverified_context
import sys, urllib, json, csv
from StringIO import StringIO

class WebServiceResource(GlycanResource):

    """
    Class for accessing glycan data from webservices

    url: Base URL for webservices

    """

    def __init__(self,*args,**kw):
        super(WebServiceResource,self).__init__(*args,**kw)
        self.attr(kw,'apiurl',required=True)
        self.attr(kw,'verbose',default=False)

    def queryws(req,url,method,kwargs):
        assert method == "GET"
	if len(kwargs) > 0:
            if '?' in url or '%(' in url:        
                url = url%kwargs
            else:
                url += "?" + urllib.urlencode(kwargs)
        # print >>sys.stderr, url
        h = urllib.urlopen(url)
        return h.read()

    def parseSection(self,name,keyvaluepairs):
        method = keyvaluepairs.get("method","GET")
        url = keyvaluepairs.get("url")
        if not url:
            req = keyvaluepairs.get("request",name)
            url = self.apiurl + "/" + req
        params = filter(None,map(lambda s: s.strip(),keyvaluepairs.get('params',"").split(',')))
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
		response = csv.DictReader(StringIO(response))
            return response

        self.set_method("query_"+str(name), _query)
        return [("query_"+name,params)]
