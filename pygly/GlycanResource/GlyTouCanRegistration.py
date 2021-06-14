
from __future__ import print_function

import sys, os.path, traceback
try:
    from urllib.parse import urlparse, urlencode
    from urllib.request import urlopen, Request, build_opener, HTTPSHandler, HTTPHandler
    from urllib.error import HTTPError
except ImportError:
    from urlparse import urlparse
    from urllib import urlencode
    from urllib2 import urlopen, Request, HTTPError, build_opener, HTTPSHandler, HTTPHandler
import json, base64

from .GlycanResource import GlycanResource

class GlyTouCanCredentialsNotFound(RuntimeError):
    pass

class GlyTouCanRegistration(GlycanResource):
    apiendpt = "https://glytoucan.org/"
    credfile = ".gtccred"
    user = None
    apikey = None
    opener = None

    def __init__(self,**kw):
        super(GlyTouCanRegistration,self).__init__(**kw)
        self.setup()
    
    def getcredentials(self):

        # script directory
        dir = os.path.split(sys.argv[0])[0]
        credfile = os.path.join(dir, GlyTouCanRegistration.credfile)
        if os.path.exists(credfile):
            user, apikey = open(credfile).read().split()
            return user, apikey

        # local directory
        credfile = GlyTouCanRegistration.credfile
        if os.path.exists(credfile):
            user, apikey = open(credfile).read().split()
            return user, apikey

        # home directory
        dir = os.path.expanduser("~")
        credfile = os.path.join(dir, GlyTouCanRegistration.credfile)
        if os.path.exists(credfile):
            user, apikey = open(credfile).read().split()
            return user, apikey

        raise GlyTouCanCredentialsNotFound()

    def setup(self, user=None, apikey=None):
        if user == None:
            user = self.user
            apikey = self.apikey
        if user == None:
            user, apikey = self.getcredentials()
        self.opener = build_opener(HTTPSHandler(),HTTPHandler())
	# print(('%s:%s'%(user, apikey)).encode('utf8'))
        self.basicauthhdr =  "Basic %s"%(base64.b64encode(('%s:%s'%(user, apikey)).encode('utf8')),) 
	# print(self.basicauthhdr)

    def register(self, sequence):
        params = json.dumps(dict(sequence=sequence))
        # POST request
        req = Request(self.apiendpt+'glycan/register', params)
        req.add_header('Content-Type', 'application/json')
        req.add_header('Accept', 'application/json')
        req.add_header("Authorization", self.basicauthhdr)
        try:
            self.wait()
            response = json.loads(self.opener.open(req).read())
            if response['status'].split()[0] == '202':
                return str(response['message'])
            print(response)
            return None
        except HTTPError as e:
            print(str(e), file=sys.stderr)
        except (ValueError, IOError) as e:
            traceback.print_exc()
        return None

