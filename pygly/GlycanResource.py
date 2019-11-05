
from ReferenceTable import ReferenceTable
import os, os.path, sys, time, traceback
from collections import defaultdict
from lockfile import FileLock
import cPickle as pickle
import gzip
import csv
import urllib

import warnings                                                                                                 
warnings.filterwarnings('ignore')

class GlycanResource(ReferenceTable):
    """
       Abstract base class for glycan resources whose data is
       accessed using triple store queries or HTTP requests

       delaytime: Waiting time between request batches, default 0.2 sec
       delaybatch: Size of request batches with no waiting time, default 1 (no batches)
       retries: Maximum number of retries: default 4

    """

    def __init__(self,**kw):
        self._delaybatch = kw.get('delaybatch',1)
        self._delaytime = kw.get('delaytime',0.2)
        self._retries = kw.get('retries',4)
        self._lastrequesttime = 0
        self._requestcount = 0

        self._cache = None
        self.attr(kw,"cachefile",default=None)
        
        super(GlycanResource,self).__init__(iniFile=kw.get('iniFile'))

    def writecache(self):
	if not self._cachefile:
	    return
        filelock = FileLock(self._cachefile)
        try:
            filelock.acquire()
            cache = self.readcache()
            for key in self._cache:
                if key not in cache:
                    cache[key] = dict()
                cache[key].update(self._cache[key])
            wh = gzip.open(self._cachefile, 'wb')
            pickle.dump(cache, wh, -1)
            wh.close()
            self._cache = cache
        finally:
            filelock.release()

    def readcache(self):
	try:
	    print >>sys.stderr, "Reading cache...."
            data = pickle.load(gzip.open(self._cachefile, 'rb'))
	    print >>sys.stderr, "Done."
	    return data
	except (IOError,ValueError,TypeError,AttributeError):
	    pass
        return dict()

    def wait(self,delay=None):
        elapsed = time.time() - self._lastrequesttime
        if delay != None:
            if elapsed < delay:
                time.sleep(delay-elapsed)
        elif (self._requestcount % self._delaybatch) == 0 and self._requestcount > 0:
            if elapsed < delay:
                time.sleep(self._delaytime-elapsed)
        self._lastrequesttime = time.time()
        self._requestcount += 1

    def attr(self,kw,key,default=None,required=False):
        if hasattr(self,key):
            setattr(self,"_"+key,getattr(self,key))
        elif key in kw:
            setattr(self,"_"+key,kw['key'])
        elif not required:
            setattr(self,"_"+key,default)
        else:
            raise RuntimeError("Can't find class/instance parameter %s for class %s"%(key,self.__class__.__name__))

import rdflib

class TripleStoreResource(GlycanResource):

    """
    Class for accessing glycan data from triple stores using arbitrary SPARQL queries

    endpt: SPARQL query endpoint
    defns: Default accession namespace prefix for URIs

    Config file provides SPARQL queryies and the names of the methods
    that should be created to access them...

    """

    def __init__(self,**kw):
        super(TripleStoreResource,self).__init__(**kw)
        self.attr(kw,'defns',default=None)
        self.attr(kw,'endpt',required=True)
        if self._defns:
            self._ns = rdflib.Namespace(self._defns)
        self._ts = rdflib.ConjunctiveGraph(store='SPARQLStore')
        self._ts.open(self._endpt)

    def queryts(self,sparql):
        self.wait()

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
        params = filter(None,map(str.strip,keyvaluepairs.get('params',"").split(',')))

        def _query(self,*args,**kw):
            # print >>sys.stderr, "query_%s: %s, %s"%(name,args,kw)
            kwargs = {}
            for i,param in enumerate(params):
                if param in keyvaluepairs:
                    kwargs[param] = keyvaluepairs[param]
                if i < len(args):
                    kwargs[param] = args[i]
                elif kw.get(param) != None:
                    kwargs[param] = kw[param]
                assert param in kwargs, " ".join(map(repr,[param, kwargs]))
            sparqlstr = sparql%kwargs
            response = self.queryts(sparqlstr)
            vars = map(str,response.vars)
            for row in response.bindings:
                row = tuple(map(row.get,response.vars))
                yield dict(zip(vars,map(str,row)))

        self.set_method("query_"+name, _query)
        return [("query_"+name,params)]

class GlyTouCan(TripleStoreResource):

    endpt = "http://ts.glytoucan.org/sparql"
    defns = "http://rdf.glycoinfo.org/glycan/"
    cachefile = ".gtccache_new"
    
    sequence_formats = set(["wurcs", "glycoct", "iupac_extended", "iupac_condensed"])

    crossref_resources = set(['glycosciences_de', 'pubchem', 'kegg',
                              'unicarbkb', 'glyconnect', 'glycome-db',
                              'unicarb-db', 'carbbank', 'pdb', 'cfg',
                              'bcsdb'])

    @staticmethod
    def partition(fn):
        def wrapper(self,**kw):
            # print >>sys.stderr, "partition: %s"%(kw,)
            if 'accession' not in kw:
                for i in range(0,10):
                    for row in fn(self,accession="G%0d.*"%(i,),**kw):
                        yield row
            else:
                for row in fn(self,**kw):
                    yield row
        return wrapper

    @staticmethod
    def prefetch(fn):
        def wrapper(self,**kw):
            # print >>sys.stderr, "prefetch:",kw
	    if not self._cache:
	        self._cache = self.readcache()
            kw1 = dict((k,v) for k,v in kw.items() if k != 'accession')
            key = fn.__name__+":"+":".join("%s=%s"%(k,v) for k,v in sorted(kw1.items()))
            # print >>sys.stderr, "cache key:",key
            if key not in self._cache:
                # print >>sys.stderr, "fill cache:",key
                self._cache[key] = dict()
                for row in fn(self,**kw1):
                    if row['accession'] not in self._cache[key]:
                        self._cache[key][row['accession']] = []
                    self._cache[key][row['accession']].append(row)
                self.writecache()
            if 'accession' in kw:
                for row in self._cache[key].get(kw['accession'],[]):
                    yield row
            else:
                for acc in self._cache[key]:
                    for row in self._cache[key][acc]:
                        yield row
        return wrapper

    def __init__(self,**kw):
        super(GlyTouCan,self).__init__(**kw)
        for k in self.keys():
            self.modify_method(k,self.partition)
	    if kw.get('usecache',False):
                self.modify_method(k,self.prefetch)

    def getseq(self,accession,format='wurcs'):
        assert format in self.sequence_formats
        for row in self.query_sequence(accession=accession,format=format):
            return row['sequence']
        return None

    def allseq(self,format=None):
        assert format == None or format in self.sequence_formats
        for row in self.query_sequence(format=format):
            yield row['accession'], row['format'], row['sequence']

    def getmass(self,accession):
        for row in self.query_mass(accession=accession):
            try:
                return float(row['mass'])
            except ValueError:
                pass
        return None

    def allmass(self):
        for row in self.query_mass():
            try:
                yield row['accession'],float(row['mass'])
            except ValueError:
                pass

    def getmonocount(self,accession):
        for row in self.query_monocount(accession=accession):
            try:
                return int(row['count'])
            except ValueError:
                pass
        return None

    def allmonocount(self):
        for row in self.query_monocount():
            try:
                yield row['accession'],int(row['count'])
            except ValueError:
                pass

    def getrefs(self,accession):
        refs = set()
        for row in self.query_references(accession=accession):
            try:
                refs.add(int(row['ref']))
            except ValueError:
                continue
        return sorted(refs)

    def allrefs(self):
        for row in self.query_references():
            try:
                yield row['accession'],int(row['ref'])
            except ValueError:
                continue

    def getcrossrefs(self,accession,resource=None):
        assert resource == None or resource in self.crossref_resources
        for row in self.query_crossrefs(accession=accession,resource=resource):
            if row['resource'] in self.crossref_resources:
                yield row['resource'],row['entry']

    def allcrossrefs(self,resource=None):
        assert resource == None or resource in self.crossref_resources
        for row in self.query_crossrefs(resource=resource):
            if row['resource'] in self.crossref_resources:
                yield row['accession'],row['resource'],row['entry']

    def getmotifs(self,accession):
        return [ row['motif'] for row in self.query_motifs(accession=accession) ]
            

    def allmotifaligns(self):
        for row in self.query_motifs():
            yield row['accession'],row['motif']

    def allmotifs(self):
        for row in self.query_allmotif():
            yield row['accession'],row['label'],row['redend']

    def exists(self,acession):
        for row in self.query_exists(accession=accession):
            return True
        return False

    def allaccessions(self):
        for row in self.query_exists():
            yield row['accession']

    def gethash(self,accession):
	for row in self.query_hash(accession=accession):
	    yield row['hash']

    def allhash(self):
	for row in self.query_hash():
	    yield row['accession'],row['hash']

    def gettaxa(self,accession):
	for row in self.query_taxonomy(accession=accession):
	    yield row['taxon']

    def bytaxa(self,taxon):
	for row in self.query_taxonomy(taxon=taxon):
	    yield row['accession']

    def alltaxa(self):
	for row in self.query_taxonomy():
	    yield row['accession'],row['taxon']

    def gettopo(self,accession):
	for row in self.query_topology(accession=accession):
	    yield row['topology']

    def alltopo(self):
	for row in self.query_topology():
	    yield row['accession'],row['topology']

    def getcomp(self,accession):
	for acc in set(list(self.gettopo(accession)) + [accession]):
	    for row in self.query_composition(accession=acc):
	        yield row['composition']

    def allcomp(self):
        topo = defaultdict(set)
        for s, t in self.alltopo():
            topo[t].add(s)
        seen = set()
        for row in self.query_composition():
	    t,c = row['accession'],row['composition']
            for s in topo[t]:
                if (s, c) not in seen:
                    seen.add((s, c))
                    yield s, c
                if (c, c) not in seen:
                    seen.add((c, c))
                    yield c, c

    def getbasecomp(self,accession):
	for acc in set(list(self.getcomp(accession))+[accession]):
	    for row in self.query_basecomposition(accession=acc):
	        yield row['basecomposition']

    def allbasecomp(self):
	comp = defaultdict(set)
	for s, c in self.allcomp():
	    comp[c].add(s)
	seen = set()
	for row in self.query_basecomposition():
	    c,bc = row['accession'],row['basecomposition']
	    for s in comp[c]:
		if (s,bc) not in seen:
		    seen.add((s,bc))
		    yield s, bc
		if (bc,bc) not in seen:
		    seen.add((bc,bc))
		    yield bc,bc

class UniCarbKBTS(TripleStoreResource):

    endpt = "http://130.56.249.35:40935/unicarbkb/query"
    defns = "http://rdf.unicarbkb.org/structure/"

    def __init__(self):
	iniFile = os.path.join(os.path.dirname(os.path.realpath(__file__)),"unicarbkb.ini")
	super(UniCarbKBTS,self).__init__(iniFile=iniFile)

    def alltaxa(self):
	for row in self.query_taxonomy():
	    yield row['accession'],row['taxon']

    def allgtc(self):
	for row in self.query_gtcacc():
	    yield row['accession'],row['glytoucan']

    def allpub(self):
	if False:
	    yield None

class UniCarbKBDump(object):
    dumpfileurl = "https://gitlab.com/matthew.campbell1980/Unicarb-Glygen/raw/master/data_files/unicarbkb/DATA_RELEASE/STABLE/%s.csv"
    species2taxa = {'human': '9606', 'mouse': '10090', 'rat': '10116'}

    def records(self):
	for species in ('human','mouse','rat'):
            url = self.dumpfileurl%(species,)
            for row in csv.DictReader(urllib.urlopen(url)):
		row['taxid'] = self.species2taxa[species]
		yield row

    def alltaxa(self):
	seen = set()
	for row in self.records():
	    data = (row['Id'],row['taxid'])
	    if data in seen:
		continue
	    seen.add(data)
	    yield data

    def allgtc(self):
	seen = set()
	for row in self.records():
	    if not row['Toucan']:
		continue
	    data = (row['Id'],row['Toucan'])
	    if data in seen:
		continue
	    seen.add(data)
	    yield data

    def allpub(self):
	seen = set()
	for row in self.records():
	    if not row['Pmid'] or row['Pmid'] == "0":
		continue
	    data = (row['Id'],row['Pmid'])
	    if data in seen:
		continue
	    seen.add(data)
	    yield data

def union(methodname):
    def wrapper(self):
	seen = set()
	for cls in (UniCarbKBTS,UniCarbKBDump):
	    try:
	        for row in getattr(cls,methodname)(self):
	            if row in seen:
		        continue
		    seen.add(row)
		    yield row
	    except IOError:
		pass
    return wrapper

class UniCarbKB(UniCarbKBDump,UniCarbKBTS):

    alltaxa = union("alltaxa")
    allgtc = union("allgtc")
    allpub = union("allpub")

if __name__ == "__main__":

    cls = sys.argv[1]
    resource = eval(cls+"()")
    query = sys.argv[2]
    method = getattr(resource,query)
    headers = None
    for r in method(*sys.argv[3:]):
	if isinstance(r,basestring):
	    print r
	elif isinstance(r,dict):
	    if headers == None:
		headers = sorted(r.keys())
		headers.remove('accession')
		headers = ['accession'] + headers
	    print "\t".join(map(r.get,headers))
	else:
	    print "\t".join(map(str,r))

    # gtc = GlyTouCan(usecache=True)
    # for acc,format,seq in gtc.allseq():
    #     print acc,format,seq
    # for acc,res,entry in gtc.allcrossrefs():
    #     print acc,res,entry
    # for i,(acc,mass) in enumerate(gtc.allmass()):
    #     print i,acc,mass
