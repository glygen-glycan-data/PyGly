
import os, os.path, sys, re
import ConfigParser
import mwclient 
from rdflib import Dataset, ConjunctiveGraph, Namespace, Literal, BNode, URIRef

class SMWCredentialsNotFound(RuntimeError):
    pass

class SMWSite(object):
    credfile = ".smwcred"
    protocol = 'http'
    host = 'localhost'
    port = 80
    tsport = 3030
    prefix = None
    propns = None

    def __init__(self,**kwargs):
        
	assert kwargs.get('prefix',self.prefix) 
	theprotocol = kwargs.get('protocol',self.protocol)
        thehost = kwargs.get('host',self.host) + ":" + str(kwargs.get('port',self.port))
	thepath = '/' + kwargs.get('prefix',self.prefix) + '/'
        self.site = mwclient.Site((theprotocol,thehost),path=thepath)
        self.site.login(*self.getcredentials())
	self.tsrw = Dataset(store='SPARQLUpdateStore')
	self.tsro = ConjunctiveGraph(store='SPARQLUpdateStore')
        baseurl = "%s://%s:%s/%s/"%(kwargs.get('protocol',self.protocol),kwargs.get('host',self.host),
                                    kwargs.get('tsport',self.tsport),kwargs.get('prefix',self.prefix))
        self.tsrw.open((baseurl+'query',baseurl+'update'))
	self.tsro.open(baseurl+'query')
        semanticns = "http://glycandata.glygen.org/"+ kwargs.get('prefix',self.prefix) +"/Special:URIResolver/"
        self.propns = Namespace(semanticns+"Property-3A")

    def getcredentials(self):

        credfiles = []

        # script directory
        dir = os.path.split(sys.argv[0])[0]
        credfile = os.path.join(dir,self.credfile)
        
        if os.path.exists(credfile):
            credfiles.append(credfile)

        # home directory
        dir = os.path.expanduser("~")
        credfile = os.path.join(dir,self.credfile)

        if os.path.exists(credfile):
            credfiles.append(credfile)

        # local directory
        credfile = self.credfile

        if os.path.exists(credfile):
            credfiles.append(credfile)

        credentials = ConfigParser.SafeConfigParser()
        if len(credentials.read(credfiles)) == 0:
            raise SMWCredentialsNotFound()

        if not credentials.has_section(self.prefix):
            raise SMWCredentialsNotFound()

        username = credentials.get(self.prefix,"username")
        password = credentials.get(self.prefix,"password")
        if not username or not password:
            raise SMWCredentialsNotFound()
            
        return username,password

    itercat_sparql = """
        PREFIX rdfs: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
        PREFIX rdf: <http://www.w3.org/2000/01/rdf-schema#>
        PREFIX category: <http://glycandata.glygen.org/glycomotif/Special:URIResolver/Category-3A>

        SELECT ?label
        WHERE {
          ?o rdfs:type category:%(category)s . 
          ?o rdf:label ?label
        }
    """
    def itercat(self,category):
	response = self.tsro.query(self.itercat_sparql%dict(category=category))
        for row in response.bindings:
	    yield row.get(response.vars[0]).toPython()

    def deletemany(self,category=None,regex=None,allpages=None,verbose=False):
	if category:
	  for pagename in self.itercat(category):
	    if verbose:
	      print >>sys.stderr, pagename
	    self.delete(pagename)
	elif regex:
	  regex = re.compile(regex)
	  for page in self.site.pages:
	    if regex.search(page.name):
	      if verbose:
	        print >>sys.stderr, page.name
	      page.delete()
	elif allpages:
	  for page in self.site.pages:
	    if page.name != 'Main_Page':
	      if verbose:
	        print >>sys.stderr, page.name
	      page.delete()

    def delete(self,pagename):
        page = self.site.pages[pagename]
        if page.exists:
            page.delete()
                     
    def get(self,name):
        page = self.site.pages[name]
        if not page.exists:
            return None
        # print repr(page.text())
        data = dict()
        splrest = map(lambda s: s.strip(),re.split(r'\|(\w+)=',page.text()))
        splrest[0]=splrest[0].lstrip('{').strip()
        splrest[-1]=splrest[-1].rstrip('}').strip()
        cls = splrest[0]
        for i in range(1,len(splrest),2):
            key = splrest[i]
            val = splrest[i+1]
            data[key] = val
	data['template'] = cls
	data['pagename'] = name
        return self.template2class[cls](**data)

    def update(self,newobj):
        pagename = newobj.pagename(**newobj.data)
        obj = self.get(pagename)
        if obj:
            obj.data.update(newobj.data)
            return self.put(obj)
        return self.put(newobj)

    uribylabel_sparql = """
        PREFIX rdf: <http://www.w3.org/2000/01/rdf-schema#>

        SELECT ?o
        WHERE {
          ?o rdf:label "%(label)s"
        }
    """
    def uribylabel(self,label):
	response = self.tsro.query(self.uribylabel_sparql%dict(label=label))
        for row in response.bindings:
	    return row.get(response.vars[0])

    def put(self,obj):
        pagename = obj.pagename(**obj.data)
	obj.set('pagename',pagename)
        page = self.site.pages[pagename]
        changed = False
        if not page.exists or page.text() != obj.astemplate():
            page.save(obj.astemplate())
            changed = True
        dummy = page.text(expandtemplates=True, cache=False)
        for key in obj.directsubmit:
            if obj.get(key):
                subj = self.uribylabel(pagename)
                pred = self.propns[obj.directsubmit[key]]
		val = Literal(obj.get(key))
                if self.tsrw.value(subj,pred) != val:
                    self.tsrw.set((subj,pred,val))
                    changed = True
        return changed

class SMWClass(object):
    def __init__(self,**kwargs):
        self.data = dict(kwargs.items())
        self.normalize()

    def normalize(self):
        return

    def __str__(self):
        l = [self.__class__.__name__]
        for k,v in sorted(self.data.items()):
            l.append("%s = %r"%(k,v))
        return "\n".join(l)

    def get(self,key,default=None):
        return self.data.get(key,default)

    def set(self,key,value):
        self.data[key] = value

    def has(self,key):
        return (key in self.data)

    def delete(self,key):
        if key in self.data:
            del self.data[key]

    def keys(self):
        return self.data.keys()

