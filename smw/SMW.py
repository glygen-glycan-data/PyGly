
import os, os.path, sys, re, inspect, copy
import ConfigParser
import mwclient 

from rdflib import Dataset, ConjunctiveGraph, Namespace, Literal, BNode, URIRef

class SMWUndefinedSiteParameter(RuntimeError):
    def __init__(self,prefix,key):
	super(SMWUndefinedSiteParameter,self).__init__("SMW site %s parameter %s not defined"%(prefix,key))

class SMWSite(object):

    def __init__(self,**kwargs):
	kwargs = self.getparams(kwargs)
	
	theprotocol = kwargs['protocol']
	thehost = kwargs['host']
        thesmwhost = thehost + ":" + kwargs['port']
	thepath = '/' + self.prefix + '/'
        self.site = mwclient.Site((theprotocol,thesmwhost),path=thepath)
	kwargs['smwurl'] = "%s://%s%s"%(theprotocol,thesmwhost,thepath)
        self.site.login(kwargs['username'],kwargs['password'])
	self.tsrw = Dataset(store='SPARQLUpdateStore')
	self.tsro = ConjunctiveGraph(store='SPARQLUpdateStore')
	thetsport = kwargs['tsport']
	baseurl = "%s://%s:%s/%s/"%(theprotocol,thehost,thetsport,self.prefix)
	kwargs['tsurl'] = baseurl
	self.tsro.open(baseurl+'query')
        self.tsrw.open((baseurl+'query',baseurl+'update'))
        ns = "http://glycandata.glygen.org/"+ self.name +"#"
        self.ns = Namespace(ns)
        subjns = "http://glycandata.glygen.org/"+ self.prefix +"/Special:URIResolver/"
        self.subjns = Namespace(subjns)
	self.params = kwargs

	self.template2class = {}
	for clsname,cls in inspect.getmembers(sys.modules[self.__class__.__module__], inspect.isclass):
	    if hasattr(cls,'template'):
		self.template2class[cls.template] = cls

    def __str__(self):
	l = ["SMW site %(prefix)s"%self.params]
	for k,v in sorted(self.params.items()):
	    if k == 'password':
	        l.append("  % 8s = %s"%(k,"X"*len(v)))
	    else:
	        l.append("  % 8s = %s"%(k,v))
	return "\n".join(l)

    def getenvironment(self,**kwargs):
	import sys
	if kwargs.get("smwenv") == "PROD":
	    return ""
	if kwargs.get("smwenv") == "DEV":
	    return "dev"
	if kwargs.get("smwenv") == "TEST":
	    return "test"
	if os.environ.get("SMWENV") == "PROD":
	    return ""
	if os.environ.get("SMWENV") == "DEV":
	    return "dev"
	if os.environ.get("SMWENV") == "TEST":
	    return "test"
        if (len(sys.argv) > 2 and sys.argv[1] == "--smwenv" and sys.argv[2] == "PROD"):
	    sys.argv.pop(1) 
	    sys.argv.pop(1) 
	    return ""
        if (len(sys.argv) > 2 and sys.argv[1] == "--smwenv" and sys.argv[2] == "DEV"):
            sys.argv.pop(1)
            sys.argv.pop(1)
	    return "dev"
        if (len(sys.argv) > 2 and sys.argv[1] == "--smwenv" and sys.argv[2] == "TEST"):
            sys.argv.pop(1)
            sys.argv.pop(1)
	    return "test"
	return "dev"

    def getparams(self,params):
	self.env = self.getenvironment(**params)
	self.name = params.get("name",self._name)
	self.prefix = self.name + self.env
	params = self.getiniparams()
	params.update(params)
	params['name'] = self.name
	params['prefix'] = self.prefix
	params['env'] = self.env
	for key in ('protocol','host','port','tsport','username','password'):
	    if key not in params:
	        raise SMWUndefinedSiteParameter(self.prefix,key)
	return params

    def getiniparams(self):

	origsmwfile = ".%s.ini"%(self.name)
        smwfiles = []

        # script directory
        dir = os.path.split(sys.argv[0])[0]
        smwfile = os.path.join(dir,origsmwfile)
        smwfiles.append(smwfile)

        # home directory
        dir = os.path.expanduser("~")
        smwfile = os.path.join(dir,origsmwfile)
        smwfiles.append(smwfile)

	# local directory

        smwfiles.append(origsmwfile)

        iniFile = ConfigParser.SafeConfigParser()
        if len(iniFile.read(smwfiles)) == 0:
            return dict()

        if not iniFile.has_section(self.prefix):
            return dict()

        return dict(iniFile.items(self.prefix))

    itercat_sparql = """
        PREFIX glycomotif: <%(ns)s>

        SELECT ?id
        WHERE {
          ?o rdf:type glycomotif:%(category)s .
          ?o glycomotif:id ?id
        }
    """
    
    def itercat(self,category):
	response = self.tsro.query(self.itercat_sparql%dict(ns=self.ns,category=category))
        for row in response.bindings:
	    yield str(row.get(response.vars[0]))

    def iternamespace(self,namespace):
        for t in self.site.allpages(namespace=namespace):
            if t.exists:
                yield t

    def itertemplates(self):
        return self.iternamespace(10)

    def itercategories(self):
        return self.iternamespace(14)

    def iterproperties(self):
        return self.iternamespace(102)

    def iterforms(self):
        return self.iternamespace(106)

    def iterpages(self,regex=None,exclude_categories=None,include_categories=None):
        if exclude_categories == None and include_categories == None and regex == None:
            for p in self.iternamespace(0):
                if p.exists:
                    yield p
        else:
            if regex != None and regex.isinstance(regex,basestring):
                regex = re.compile(regex)
            if exclude_categories != None:
                exclude_categories = set(exclude_categories)
            if include_categories != None:
                include_categories = set(include_categories)
            for p in self.iternamespace(0):
                if p.exists:
                    if regex and not regex.search(p.name):
                        continue
                    cats = set(map(lambda c: str(c.name.split(':',1)[1]),p.categories()))
                    if exclude_categories != None and len(cats&exclude_categories) > 0:
                        continue
                    if include_categories != None and len(cats&include_categories) == 0:
                        continue
                    yield p

    def dumpiterable(self,dir,iter):
        try:
            os.makedirs(dir)
        except OSError:
            pass
        for p in iter:
            name = p.name
            if ":" in name:
                name = name.split(":",1)[1]
            try:
                print p.name
            except:
                continue
            try:
                h = open(dir+'/'+name+'.txt','w')
                h.write(p.text())
            except:
                pass
            finally:
                h.close()

    def dumpsite(self,dir):
        try:
            os.makedirs(dir)
        except OSError:
            pass
        self.dumpiterable(os.path.join(dir,'templates'),self.itertemplates())
        self.dumpiterable(os.path.join(dir,'categories'),self.itercategories())
        self.dumpiterable(os.path.join(dir,'forms'),self.iterforms())
        self.dumpiterable(os.path.join(dir,'properties'),self.iterproperties())
        self.dumpiterable(os.path.join(dir,'pages'),
                          self.iterpages(exclude_categories=self.dump_exclude_categories))

    def loadsite(self,dir):
        for ns,subdir in (('MediaWiki','mediawiki'),
                          ('Property','properties'),
                          ('Category','categories'),
                          ('Template','templates'), 
                          ('Form','forms'),
                          ('','pages')):
	  if ns:
              ns += ":"
	  for f in os.listdir(os.path.join(dir,subdir)):
              if not f.endswith('.txt'):
                  continue
              pagetext = open(os.path.join(dir,subdir,f)).read()
              pagename = ns+f[:-4]
              page = self.site.pages[pagename]
	      
              if page.text() != pagetext:
                  page.save(pagetext)
                  print >>sys.stderr, pagename
        
    def deletemany(self,category=None,regex=None,allpages=None,verbose=False):
	if category:
	  for pagename in self.itercat(category):
	    if verbose:
	      print >>sys.stderr, pagename
	    self.delete(pagename)
          for p in self.iterpages(include_categories=[category]):
            if verbose:
              print >>sys.stderr, p.name
            p.delete()
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
      subj = self.uribyid(pagename)
      if subj:
        for t in self.tsrw.triples((subj,None,None)):
          self.tsrw.remove(t)
                     
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
        data['id'] = name
        return self.template2class[cls](**data)

    def update(self,newobj):
        pagename = newobj.pagename(**newobj.data)
        obj = self.get(pagename)
        if obj:
            obj.data.update(newobj.data)
            return self.put(obj)
        return self.put(newobj)

    uribyid_sparql = """
        PREFIX glycomotif: <%(ns)s>
        
        SELECT ?o
        WHERE {
          ?o glycomotif:id "%(id)s"
        }
        """
    
    def uribyid(self,id):
      response = self.tsro.query(self.uribyid_sparql%dict(ns=self.ns,id=id))
      for row in response.bindings:
        return row.get(response.vars[0])
      # Hack - the URI isn't in the database yet, so we have to scramble.
      # Essentially we need to simulate the MW label to URI algorithm
      if re.search(r'^[A-Za-z0-9.]+$',id):
        return self.subjns[id]
      return None

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
                subj = self.uribyid(pagename)
                pred = self.ns[obj.directsubmit[key]]
		val = Literal(obj.get(key))
                if self.tsrw.value(subj,pred) != val:
                    self.tsrw.set((subj,pred,val))
                    changed = True
        return changed

class SMWClass(object):
    directsubmit = []
    def __init__(self,**kwargs):
        self.data = self.toPython(kwargs)

    @staticmethod
    def pagename(**kwargs):
        assert kwargs.get('id')
        return kwargs.get('id')

    @staticmethod
    def asboolean(value):
	if isinstance(value,basestring):
	    if value.lower() in ("true","1","yes","t","y"):
		return True
	    elif value.lower() in ("false","0","no","f","n"):
		return False
	    else:
		raise ValueError("Can't infer boolean value from %r"%(value,))
	if not isinstance(value,bool):
	    raise TypeError("Can't infer boolean value from %r"%(value,))
	return value

    def toPython(self,data):
        for k in list(data.keys()):
            if data.get(k) in (None,"",[]):
                del data[k]
	return data

    def toTemplate(self,data):
	for k in list(data.keys()):
	    if k in ('template','id','pagename'):
		del data[k]
		continue
	    if data.get(k) in (None,"",[]):
		del data[k]
		continue
	return data

    def astemplate(self):
        l = [ "{{%s"%(self.template,) ] 
	data = self.toTemplate(copy.copy(self.data))
        for k in sorted(data.keys()):
            v = data.get(k)
            l.append("|%s=%s"%(k,v))
        l.append("}}")
        return "\n".join(l)

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

