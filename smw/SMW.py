
import os, os.path, sys, re, inspect, copy, traceback
import ConfigParser
import mwclient 

from rdflib import Dataset, ConjunctiveGraph, Namespace, Literal, BNode, URIRef

import requests

class SMWUndefinedSiteParameter(RuntimeError):
    def __init__(self,prefix,key):
	super(SMWUndefinedSiteParameter,self).__init__("SMW site %s parameter %s not defined"%(prefix,key))

class SMWSite(object):

    def __init__(self,**kwargs):
	defaults = dict(protocol='http',port=None)
	kwargs = self.getparams(kwargs,defaults)
	theprotocol = kwargs['protocol']
	thesmwhost = kwargs['host']
	if kwargs.get('port') != None:
            thesmwhost += (":" + kwargs['port'])
	thepath = '/'
        try:
	    self.site = mwclient.Site(scheme=theprotocol,host=thesmwhost,path=thepath,reqs=dict(timeout=300))
	except TypeError:
	    #mwclient, pre 0.10.0
            self.site = mwclient.Site(host=(theprotocol,thesmwhost),path=thepath)
	kwargs['smwurl'] = "%s://%s%s"%(theprotocol,thesmwhost,thepath)
        self.site.login(kwargs['username'],kwargs['password'])
	self.params = kwargs

	self.template2class = {}
	for clsname,cls in inspect.getmembers(sys.modules[self.__class__.__module__], inspect.isclass):
            if 'template' in cls.__dict__:
		self.template2class[cls.template] = cls
        if not kwargs.get('quiet',False):
	    print >>sys.stderr, "Connected to %s"%(self.title())

    def title(self):
	env = ""
        if self.env:
	    env = " (%s)"%(self.env.upper(),)
	return self.name+env

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
	if os.environ.get("SMWENV") == "PROD":
	    return ""
	if os.environ.get("SMWENV") == "DEV":
	    return "dev"
	if os.environ.get("SMWENV") == "TEST":
	    return "test"
	return "dev"

    def getparams(self,kwargs,defaults={}):
	self.env = self.getenvironment(**kwargs)
	self.name = kwargs.get("name",self._name)
	self.prefix = self.name + self.env
	params = self.getiniparams()
	params.update(kwargs)
	params['name'] = self.name
	params['prefix'] = self.prefix
	params['env'] = self.env
	for key in ('protocol','host','port','username','password'):
	    if key not in params:
		if key not in defaults:
	            raise SMWUndefinedSiteParameter(self.prefix,key)
		else:
		    params[key] = defaults[key]
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

    def itercat(self,category):
	for p in self.iterpages(include_categories=[category]):
	    yield p.name

    def iterall(self):
	for p in self.itermediawiki():
	    yield p
	for p in self.itertemplates():
	    yield p
	for p in self.itercategories():
	    yield p
	for p in self.iterforms():
	    yield p
	for p in self.iterpages():
	    yield p
	for p in self.itertalk():
	    yield p

    def iternamespace(self,namespace):
        for t in self.site.allpages(namespace=namespace):
            if t.exists:
                yield t

    def itermediawiki(self):
        for p in self.iternamespace(8):
	    if p.name in ('MediaWiki:Smw import foaf','MediaWiki:Smw import owl','MediaWiki:Smw import skos'):
                continue
            yield p

    def itertemplates(self):
        return self.iternamespace(10)

    def itercategories(self):
        return self.iternamespace(14)

    def iterproperties(self):
        for p in self.iternamespace(102):
	    if p.name.startswith('Property:Foaf:') or p.name.startswith('Property:Owl:'):
		continue
	    yield p

    def iterforms(self):
        return self.iternamespace(106)

    def itertalk(self):
        return self.iternamespace(1)

    def iterpages(self,regex=None,exclude_regex=None,exclude_categories=None,include_categories=None):
        if exclude_categories == None and include_categories == None and regex == None and exclude_regex == None:
            for p in self.iternamespace(0):
                if p.exists:
                    yield p
        else:
            if regex != None and isinstance(regex,basestring):
                regex = re.compile(regex)
            if exclude_regex != None and isinstance(exclude_regex,basestring):
                exclude_regex = re.compile(exclude_regex)
            if exclude_categories != None:
                exclude_categories = set(exclude_categories)
            if include_categories != None:
                include_categories = set(include_categories)
            for p in self.iternamespace(0):
                if p.exists:
                    if regex and not regex.search(p.name):
                        continue
		    if exclude_regex and exclude_regex.search(p.name):
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
	    if not dir.endswith('/pages'):
                name = name.split(":",1)[1]
            try:
                print p.name
            except:
                continue
            try:
                h = open(dir+'/'+name+'.txt','w')
                h.write((p.text().rstrip('\n')+'\n').encode('utf8').replace(r'\n',r'\\n').decode('string_escape'))
            except:
                raise
            finally:
                h.close()

    def dumpsite(self,dir,exclude_regex=None,exclude_categories=set()):
        try:
            os.makedirs(dir)
        except OSError:
            pass
        self.dumpiterable(os.path.join(dir,'mediawiki'),self.itermediawiki())
        self.dumpiterable(os.path.join(dir,'templates'),self.itertemplates())
        self.dumpiterable(os.path.join(dir,'categories'),self.itercategories())
        self.dumpiterable(os.path.join(dir,'forms'),self.iterforms())
        self.dumpiterable(os.path.join(dir,'properties'),self.iterproperties())
        self.dumpiterable(os.path.join(dir,'pages'),
                          self.iterpages(exclude_regex=exclude_regex,
                                         exclude_categories=exclude_categories))
        self.dumpiterable(os.path.join(dir,'talk'),self.itertalk())

    def loadsite(self,dir):
        for ns,subdir in (('MediaWiki','mediawiki'),
                          ('Property','properties'),
                          ('Category','categories'),
                          ('Template','templates'), 
                          ('Form','forms'),
                          ('','pages'),
			  ('Talk','talk')):
	  if not os.path.exists(os.path.join(dir,subdir)):
	      continue
	  if ns:
              ns += ":"
	  for f in os.listdir(os.path.join(dir,subdir)):
              if not f.endswith('.txt'):
                  continue
              pagetext = open(os.path.join(dir,subdir,f)).read().rstrip().decode('utf8')
              pagename = ns+f[:-4]
              page = self.site.pages[pagename]
	      
              if page.text() != pagetext:
                  print >>sys.stderr, pagename
                  page.save(pagetext)
        
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

    # def update(self,**kw):
    #   self.data.update(kw)
 
    def delete(self,pagename):
      page = self.site.pages[pagename]
      if page.exists:
        page.delete()

    def _get(self,name):
	return self.site.pages[name]

    @staticmethod
    def parse_template_text(name,thepage):
	obj = None
        chunks = []
        level = 0
	lastpos = 0
        for m in re.finditer(r'({{)|(}})', thepage):
            if m.group() == '{{':
                if level == 0:
                    startpos = m.end()
		# print "  "*level,thepage[lastpos:m.end()]
		lastpos = m.end()
                level += 1
            if m.group() == '}}':
		# print "  "*level,thepage[lastpos:m.start()]
		lastpos = m.start()
                level -= 1
                if level == 0:
                    endpos = m.start()
		    # print "!",thepage[startpos:endpos]
                    chunks.append(thepage[startpos:endpos])
        assert level == 0, name

        subobjs = []
        for index, chunk in enumerate(chunks[1:]+chunks[:1]):
            data = dict()
            splrest = map(lambda s: s.strip(),re.split(r'\|(\w+)=',chunk))
            cls = splrest[0]
            for i in range(1,len(splrest),2):
                key = splrest[i]
                val = splrest[i+1]
                data[key] = val
            data['template'] = cls
            data['pagename'] = name
            data['id'] = name
            if index == len(chunks)-1:
                if len(subobjs) > 0:
                    data['_subobjs'] = subobjs
                obj = data
            else:
                subobjs.append(data)
            
        return obj
        
                     
    def get(self,name):
	name = name.strip()
	if not name:
	    return None
        page = self._get(name)
        if not page.exists:
            return None

        thepage = page.text()
        obj = self.parse_template_text(name,thepage)

        if '_subobjs' in obj:
            for i in range(len(obj['_subobjs'])):
                subobj = obj['_subobjs'][i]
                cls = subobj['template']
                obj['_subobjs'][i] = self.template2class[cls](**subobj)
        cls = obj['template']
        obj = self.template2class[cls](**obj)
        
        return obj

    def update(self,newobj):
        pagename = newobj.pagename(**newobj.data)
        obj = self.get(pagename)
        if obj:
            thedata = dict(obj.data,**newobj.data)
            if thedata == obj.data:
                return False
	    obj.data = thedata
            return self.put(obj)
        return self.put(newobj)

    def put(self,obj):
        pagename = obj.pagename(**obj.data)
	obj.set('pagename',pagename)
        page = self.site.pages[pagename]
        changed = False
        if not page.exists or page.text() != obj.astemplate():
            page.save(obj.astemplate())
            changed = True
        return changed

    def refresh(self,page):
	if isinstance(page,SMWClass):
	    pagename = page.pagename(**page.data)
	    page = self.site.pages[pagename]
	elif isinstance(page,basestring):
	    page = self.site.pages[page]
	if page.exists:
	    page.purge()
            dummy = self.site.get('parse',page=page.name,prop="text")

class SMWClass(object):

    directsubmit = []
    def __init__(self,**kwargs):
	if kwargs.get('template'):
	    self.template = kwargs.get('template')
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

    @staticmethod
    def smwescape(value):
        if not re.search(r'[]{}[=|<]',value):
            return value.strip()
        # value = value.replace('{','{{(}}')
        # value = value.replace('}','{{)}}')
	value = re.sub(r'[{}]', lambda mo: '{{(}}' if mo.group() == '{' else '{{)}}', value)
        value = value.replace('|','{{!}}')
        value = value.replace('=','{{=}}')
        value = value.replace('[','{{!(}}')
        value = value.replace(']','{{!)}}')
        value = value.replace('<','{{lessthan}}')
        return value.strip()

    @staticmethod
    def smwunescape(value):
        if not re.search(r'\{\{.*\}\}',value):
            return value.strip()
        value = value.replace('{{!}}','|')
        value = value.replace('{{=}}','=')
        value = value.replace('{{(}}','{')
        value = value.replace('{{)}}','}')
        value = value.replace('{{!(}}','[')
        value = value.replace('{{!)}}',']')
        value = value.replace('{{lessthan}}','<')
        return value.strip()

    def toPython(self,data):
        for k in list(data.keys()):
            if data.get(k) in (None,[]):
                del data[k]
                continue
            if isinstance(data[k],basestring):
                data[k] = self.smwunescape(data[k])
		if not data[k]:
		    del data[k]
		    continue
	return data

    def toTemplate(self,data):
	for k in list(data.keys()):
	    if k in ('template','id','pagename'):
		del data[k]
		continue
	    if data.get(k) in (None,"",[]):
		del data[k]
		continue
            if isinstance(data[k],basestring):
                data[k] = self.smwescape(data[k])
	return data

    def astemplate(self):
        l = [ "{{%s"%(self.template,) ] 
	data = self.toTemplate(copy.copy(self.data))
        for k in sorted(data.keys()):
            if k != "_subobjs":
                v = data.get(k)
                l.append("|%s=%s"%(k,v))
        l.append("}}")
        if "_subobjs" in data:
            for subobj in data['_subobjs']:
                l[-1] += subobj.astemplate()
        return "\n".join(l)

    def __str__(self):
        l = [self.__class__.__name__]
        for k,v in sorted(self.data.items()):
            if isinstance(v,list) and isinstance(v[0],SMWClass):
                l.append("%s = ["%k)
                for i,vi in enumerate(v):
                    si = str(vi).strip()
                    si = "  " + si.replace('\n','\n  ')
                    if i > 0:
                        si = "\n" + si;
                    if i < (len(v)-1):
                        si += ','
                    l.append(si)
                l.append("  ]")
            else:
                l.append("%s = %r"%(k,v))
        return "\n".join(l)

    def save(self,wiki):
        return wiki.put(self)

    def get(self,key,default=None):
        return self.data.get(key,default)

    def set(self,key,value):
        self.data[key] = value

    def has(self,key):
        return (key in self.data)

    def append(self,key,value):
        if key not in self.data:
            self.data[key] = []
        self.data[key].append(value)

    def add(self,key,value):
        if key not in self.data:
            self.data[key] = set()
        self.data[key].add(value)

    def sort(self,key,keyfn):
	if key in self.data:
	    self.data[key].sort(key=keyfn)

    def update(self,**kwargs):
        self.data.update(self.toPython(kwargs))

    def delete(self,key):
        if key in self.data:
            del self.data[key]

    def keys(self):
        return self.data.keys()

    def items(self):
        return self.data.items()
