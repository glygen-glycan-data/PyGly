
import sys
import re
import traceback

import warnings                                                                                                              
warnings.filterwarnings('ignore')

import rdflib

import shelve

try:
    from pygly.lockfile import FileLock
except ImportError:
    from lockfile import FileLock

from GlycanResource import GlycanResource

class TripleStoreResource(GlycanResource):

    """
    Class for accessing glycan data from triple stores using arbitrary SPARQL queries

    endpt: SPARQL query endpoint
    defns: Default accession namespace prefix for URIs

    Config file provides SPARQL queryies and the names of the methods
    that should be created to access them...

    """

    def __init__(self,*args,**kw):
        super(TripleStoreResource,self).__init__(*args,**kw)
        self.attr(kw,'defns',default=None)
        self.attr(kw,'endpt',required=True)
	self.attr(kw,'verbose',default=False)

        self.attr(kw,'cachefile',default=False)
	self.attr(kw,'usecache',default=False)
        self.attr(kw,'cachemode',default='r')
	self.attr(kw,'prefetch',default=False)
        if not self._prefetch:
            self._usecache = False

        if self._defns:
            self._ns = rdflib.Namespace(self._defns)
        self._ts = rdflib.ConjunctiveGraph(store='SPARQLStore')
	self._ts.store.method = 'POST'
        self._ts.open(self._endpt)

        self._cache = None

        if self._prefetch:
            self._cache = {}

        if self._usecache:
            self.opencache()

        if self._verbose:
	    print >>sys.stderr, "TripleStoreResource:prefetch = %s"%(self._prefetch)
	    print >>sys.stderr, "TripleStoreResource:usecache = %s"%(self._usecache)
	    print >>sys.stderr, "TripleStoreResource:cachemode = %s"%(self._cachemode)

    def __del__(self):
        if hasattr(self,'_usecache') and self._usecache:
            self.closecache()

    def opencache(self):
	if self._cachefile:
            if self._verbose:
                print >>sys.stderr, "Opening cachefile:",self._cachefile
	    self._cacheondisk = shelve.open(self._cachefile,flag=self._cachemode)
	    self._cachedirty = {}
        
    def writecache(self):
	if not self._cachefile:
	    return
        filelock = FileLock(self._cachefile)
        try:
            filelock.acquire()
	    for key in self._cache:
		if self._cachedirty.get(key,False):
                    if self._verbose:
                        print >>sys.stderr, "Write key",key,"to cachefile:",self._cachefile
		    # write out "query" for key
		    self._cacheondisk[key] = self._cache[key]
		    self._cachedirty[key] = False
	    self._cacheondisk.sync()
        finally:
            filelock.release()

    def closecache(self):
        if hasattr(self,'_verbose') and hasattr(self,'_cachefile') and self._verbose:
            print >>sys.stderr, "Closing cachefile:",self._cachefile
        if hasattr(self,'_cacheondisk'):
            self._cacheondisk.close()

    def queryts(self,sparql):
        self.wait()

	if self._verbose:
	    print >>sys.stderr, "SPARQL Query:\n"
	    print >>sys.stderr, sparql
	    print >>sys.stderr, ""

        attempt = 0
        response = None
        delay = self._delaytime
        while response == None and (attempt-1) < self._retries:
            try:
                attempt += 1
                response = self._ts.query(sparql)
            except:
                traceback.print_exc()
                delay *= 2
                self.wait(delay)

        if response == None:
            raise IOError("Cannot query SPARQL endpoint")

        return response

    def triples(self,acc):
        self.wait()
        
        if not acc.startswith('http'):
            assert self._ns != None
            uri = self._ns[acc]
        else:
            uri = rdflib.Namespace(acc)[""]

        seen = set()
        for subj,pred,obj in self._ts.triples((uri,None,None)):
            if (subj,pred,obj) in seen:
                continue
            seen.add((subj,pred,obj))
            yield tuple(map(str,(subj,pred,obj)))
        for subj,pred,obj in self._ts.triples((None,None,uri)):
            if (subj,pred,obj) in seen:
                continue
            seen.add((subj,pred,obj))
            yield tuple(map(str,(subj,pred,obj)))

    def set_method(self,name,func):
        setattr(self.__class__, name, func)
        func.__name__ = name
    
    def modify_method(self,name,func):
        newfunc = func(getattr(self.__class__,name))
        setattr(self.__class__, name, newfunc)
        newfunc.__name__ = name
    
    def parseSection(self,name,keyvaluepairs):
        sparql = keyvaluepairs['sparql']
        params = filter(None,map(lambda s: s.strip(),keyvaluepairs.get('params',"").split(',')))
  	escape = filter(None,map(lambda s: s.strip(),keyvaluepairs.get('escape',"").split(',')))

        def _query(self,*args,**kw):
            # print >>sys.stderr, "query_%s: %s, %s"%(name,args,kw)
            kwargs = {}
            for i,param in enumerate(params):
                if param in keyvaluepairs:
                    kwargs[param] = keyvaluepairs[param]
                if i < len(args):
                    kwargs[param] = args[i]
		    if param in escape:
			kwargs[param] = re.escape(kwargs[param]).replace("\\","\\\\")
                elif kw.get(param) != None:
                    kwargs[param] = kw[param]
		    if param in escape:
			kwargs[param] = re.escape(kwargs[param]).replace("\\","\\\\")
                assert param in kwargs, " ".join(map(repr,[param, kwargs]))
            sparqlstr = sparql%kwargs
            response = self.queryts(sparqlstr)
            vars = map(str,response.vars)
            for row in response.bindings:
                row = tuple(map(row.get,response.vars))
                yield dict(zip(vars,map(self.tostr,row)))

        self.set_method("query_"+name, _query)
        return [("query_"+name,params)]

    @staticmethod
    def tostr(value):
        if value == None:
            return None
	try:
            return str(value)
	except ValueError:
	    pass
        return value
