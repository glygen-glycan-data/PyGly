
import sys, os.path
import urllib2
import json

from GlycanResource import GlycanResource

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
        # print user,apikey
        self.password_mgr = urllib2.HTTPPasswordMgrWithDefaultRealm()
        self.password_mgr.add_password(None, self.apiendpt, user, apikey)
        self.handler = urllib2.HTTPBasicAuthHandler(self.password_mgr)
        # self.opener = urllib2.build_opener(self.handler,urllib2.HTTPSHandler(debuglevel=1),urllib2.HTTPHandler(debuglevel=1))
        self.opener = urllib2.build_opener(self.handler)

    def register(self, sequence):
        params = json.dumps(dict(sequence=sequence))
	# POST request
        req = urllib2.Request(self.apiendpt+'glycan/register', params)
        req.add_header('Content-Type', 'application/json')
        req.add_header('Accept', 'application/json')
        try:
            self.wait()
            response = json.loads(self.opener.open(req).read())
            return response
        except (ValueError, IOError), e:
	    pass
        return None

