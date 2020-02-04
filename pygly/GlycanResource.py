
from ReferenceTable import ReferenceTable
import os, os.path, sys, time, traceback, re
from dateutil import parser as dateutil_parser
from collections import defaultdict
from lockfile import FileLock
import shelve
import cPickle as pickle
import gzip
import csv
import urllib
import json

from GlycanFormatter import WURCS20Format, GlycoCTFormat, GlycanParseError, ZeroPlusLinkCountError, UndeterminedLinkCountError, CircularError, LinkCountError
from WURCS20MonoFormatter import WURCS20MonoFormat, UnsupportedSkeletonCodeError, UnsupportedSubstituentError, InvalidMonoError

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

        self.attr(kw,"cachefile",default=None)
	if self._cachefile:
	    self._cacheondisk = shelve.open(self._cachefile,flag='c')
	    self._cache = {}
	    self._cachedirty = {}
        
        super(GlycanResource,self).__init__(iniFile=kw.get('iniFile'))

    def writecache(self):
	if not self._cachefile:
	    return
        filelock = FileLock(self._cachefile)
        try:
            filelock.acquire()
	    for key in self._cache:
		if self._cachedirty.get(key,False):
		    # write out "query" for key
		    self._cacheondisk[key] = self._cache[key]
		    self._cachedirty[key] = False
        finally:
            filelock.release()

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
            setattr(self,"_"+key,kw[key])
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

    def __init__(self,*args,**kw):
        super(TripleStoreResource,self).__init__(*args,**kw)
        self.attr(kw,'defns',default=None)
        self.attr(kw,'endpt',required=True)
	self.attr(kw,'verbose',default=False)
        if self._defns:
            self._ns = rdflib.Namespace(self._defns)
        self._ts = rdflib.ConjunctiveGraph(store='SPARQLStore')
        self._ts.open(self._endpt)

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

# Defaults ensure a 10-way partition of GlyTouCan accessions
def partitioner(kwarg="accession",fmt="G%%0%dd.*",digits=1):
    fmtstr = fmt%(digits,)
    def partition(fn):
        def wrapper(self,*args,**kw):
            if kwarg not in kw:
                for i in range(0,10**digits):
                    for row in fn(self,*args,accession=fmtstr%(i,),**kw):
                        yield row
            else:
                for row in fn(self,*args,**kw):
                    yield row
        return wrapper
    return partition

class GlyTouCanTS(TripleStoreResource):

    endpt = "http://ts.glytoucan.org/sparql"
    defns = "http://rdf.glycoinfo.org/glycan/"
    cachefile = ".gtccache_new"
    # verbose = True
    
    sequence_formats = set(["wurcs", "glycoct", "iupac_extended", "iupac_condensed"])

    crossref_resources = set(['glycosciences_de', 'pubchem', 'kegg',
                              'unicarbkb', 'glyconnect', 'glycome-db',
                              'unicarb-db', 'carbbank', 'pdb', 'cfg',
                              'bcsdb','matrixdb','glycoepitope'])

    def __init__(self,*args,**kwargs):
	kwargs['iniFile'] = os.path.join(os.path.dirname(os.path.realpath(__file__)),"glytoucan.ini")
	super(GlyTouCanTS,self).__init__(*args,**kwargs)

    @staticmethod
    def prefetch(fn):
        def wrapper(self,**kw):
            kw1 = dict((k,v) for k,v in kw.items() if k != 'accession')
            key = fn.__name__+":"+":".join("%s=%s"%(k,v) for k,v in sorted(kw1.items()))
            # print >>sys.stderr, "cache key:",key
            if key not in self._cache:
		if not self._cacheondisk.has_key(key):
                    # print >>sys.stderr, "fill cache:",key
		    self._cache[key] = {}
		    self._cachedirty[key] = True
                    for row in fn(self,**kw1):
                        if row['accession'] not in self._cache[key]:
                            self._cache[key][row['accession']] = []
                        self._cache[key][row['accession']].append(row)
                    self.writecache()
		else:
                    self._cache[key] = self._cacheondisk[key]
		    self._cachedirty[key] = False
            if 'accession' in kw:
                for row in self._cache[key].get(kw['accession'],[]):
                    yield row
            else:
                for acc in self._cache[key]:
                    for row in self._cache[key][acc]:
                        yield row
        return wrapper

    def __init__(self,**kw):
        super(GlyTouCanTS,self).__init__(**kw)
        for k in self.keys():
            self.modify_method(k,partitioner())
	    if kw.get('usecache',True):
                self.modify_method(k,self.prefetch)

    def getseq(self,accession,format='wurcs'):
        assert format in self.sequence_formats
        for row in self.query_sequence(accession=accession,format=format):
	    if row['format'] == 'wurcs' and not row['sequence'].startswith('WURCS'):
		continue
	    if row['format'] == 'glycoct' and not row['sequence'].startswith('RES'):
		continue
            return row['sequence']
        return None

    def allseq(self,format=None):
        assert format == None or format in self.sequence_formats
        for row in self.query_sequence(format=format):
	    if row['format'] == 'wurcs' and not row['sequence'].startswith('WURCS'):
                continue
	    if row['format'] == 'glycoct' and not row['sequence'].startswith('RES'):
                continue
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

    # Named for consistency with GlyTouCan class...
    def getmotif(self,accession):
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

    def allinvalid(self):
	for row in self.query_invalid():
	    yield row['accession']

    def invalid(self,accession):
	for row in self.query_invalid(accession=accession):
	    return True
	return False

    def allaccessions(self):
        for row in self.query_exists():
            yield row['accession']

    def gethash(self,accession):
	for row in self.query_hash(accession=accession):
	    return row['hash']

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
	    return row['topology']
	return None

    def alltopo(self):
	for row in self.query_topology():
	    yield row['accession'],row['topology']

    def getcomp(self,accession):
	for acc in set([self.gettopo(accession),accession]):
	    for row in self.query_composition(accession=acc):
	        return row['composition']

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
	for acc in set([self.getcomp(accession),accession]):
	    for row in self.query_basecomposition(accession=acc):
	        return row['basecomposition']

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

    def _query_date_helper(self,**kwargs):
        lastacc = None; lastaccdates = set()
	for row in sorted(self.query_date(**kwargs),key=lambda r: r['accession']):
	    row['date'] = dateutil_parser.parse(row['date']).date()
	    if row['accession'] != lastacc:
		if len(lastaccdates) > 0:
		    yield lastacc,min(lastaccdates).isoformat(),max(lastaccdates).isoformat()
		lastacc = row['accession']
		lastaccdates = set([row['date']])
	    else:
		lastaccdates.add(row['date'])
	if len(lastaccdates) > 0:
	    yield lastacc,min(lastaccdates).isoformat(),max(lastaccdates).isoformat()

    def getdate(self,accession):
	for row in self._query_date_helper(accession=accession):
	    yield row[1],row[2]

    def alldate(self):
	for row in self._query_date_helper():
	    yield row

class GlyTouCanUtil(object):
    _wurcs_mono_format = WURCS20MonoFormat()
    _wurcs_format = WURCS20Format()
    _glycoct_format = GlycoCTFormat()

    def getUnsupportedCodes(self, acc):
        codes = set()
	substs = set()
	invalid = set()
	other = set()
        sequence = self.getseq(acc, 'wurcs')
        if not sequence:
            return codes, substs, invalid, other
        monos = sequence.split('/[', 1)[1].split(']/')[0].split('][')
        for m in monos:
            try:
                g = self._wurcs_mono_format.parsing(m)
            except UnsupportedSkeletonCodeError, e:
                codes.add(e.message.rsplit(None, 1)[-1])
            except UnsupportedSubstituentError, e:
                substs.add(e.message.rsplit(None, 1)[-1])
            except InvalidMonoError, e:
                invalid.add(e.message.rsplit(None, 1)[-1])
            except GlycanParseError:
                pass
	try:
	    g = self._wurcs_format.toGlycan(sequence)
	except ZeroPlusLinkCountError:
	    other.add("0+ link count")
	except UndeterminedLinkCountError:
	    other.add("undetermined link count")
	except CircularError:
	    other.add("circular")
	except LinkCountError:
	    other.add("bad link count")
	except GlycanParseError:
	    pass
        return codes, substs, invalid, other

    def getGlycan(self, acc, format=None):
        if not format or (format == 'wurcs'):
            sequence = self.getseq(acc, 'wurcs')
            if sequence:
                try:
                    return self._wurcs_format.toGlycan(sequence)
                except GlycanParseError:
                    pass  # traceback.print_exc()
        if not format or (format == 'glycoct'):
            sequence = self.getseq(acc, 'glycoct')
            if sequence:
                try:
                    return self._glycoct_format.toGlycan(sequence)
                except GlycanParseError:
                    pass
        return None

    def glycoct(self, acc, fetch=None):
	g = self.getGlycan(acc,fetch)
	if not g:
	    return None
	return g.glycoct()

    def umw(self, acc, fetch=None):
	g = self.getGlycan(acc,fetch)                                                                             
	if not g:
	    return None
        return g.underivitized_molecular_weight()

    def wurcs2glycoct(self, acc):
	sequence = self.getseq(acc,'wurcs')
	if sequence:
	    sequence1 = urllib.quote_plus(sequence)
	    url = 'https://api.glycosmos.org/glycanformatconverter/2.3.2-snapshot/wurcs2glycoct/'+sequence1
	    try:
	        data = json.loads(urllib.urlopen(url).read())
	        if 'GlycoCT' in data:
	            return data['GlycoCT']
	    except ValueError:
		pass
	return None

    def subsumptionbyapi(self, acc):
	sequence = self.getseq(acc,'wurcs')
	if sequence:
	    sequence1 = urllib.quote_plus(sequence)
	    url = 'https://api.glycosmos.org/subsumption/0.2.0/'+sequence1
	    data = urllib.urlopen(url).read()
	    seen = set()
	    lasts = None
	    for triple in sorted(map(lambda t: tuple(map(str.strip,map(str,map(t.get,("S","P","O"))))),json.loads(data))):
		if triple in seen:
		    continue
		seen.add(triple)
		if triple[0] != lasts:
		    if lasts != None:
			print ""
		    print triple[0]
		    lasts = triple[0]
		if triple[2] == sequence:
		    print ">>  "+"\t".join(triple[1:])
		else:
		    print "    "+"\t".join(triple[1:])

    def findskel(self, skel, maxcount=None):
	if maxcount != None:
	    maxcount = int(maxcount)
        
	for acc, format, wurcs in self.allseq(format='wurcs'):
            glycoct = self.getseq(acc,format='glycoct')
            if not glycoct:
                continue
	    monos = wurcs.split('/[', 1)[1].split(']/')[0].split('][')
	    if maxcount != None and len(monos) > maxcount:
		continue
            for mono in monos:
		msk = re.search(r'^(.*?)([-_].*)?$',mono).group(1)
		assert msk
		m = re.search(r"^%s$"%(skel,),msk)
		if m:
		    yield acc, m.group(0)

    def multiseq(self):
	counts = defaultdict(set)
	for acc,fmt,seq in self.allseq():
	    counts[(acc,fmt)].add(seq)
	for k,v in counts.items():
	    if len(v) > 1:
		yield k

class GlyTouCan(GlyTouCanTS,GlyTouCanUtil):
    pass

class GlyTouCanNoCache(GlyTouCan):
    def __init__(self):
	iniFile = os.path.join(os.path.dirname(os.path.realpath(__file__)),"glytoucan.ini")
	super(GlyTouCanNoCache,self).__init__(usecache=False,iniFile=iniFile)

class UniCarbKBTS(TripleStoreResource):

    # endpt = "http://130.56.249.35:40935/unicarbkb/query"
    endpt = "http://203.101.226.128:40935/unicarbkb/query"
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
    species2filename = {'human': 'human17122019', 'mouse': 'mouse03122019', 'rat': 'rat03122019'}

    def records(self):
	for species in ('human','mouse','rat'):
            url = self.dumpfileurl%(self.species2filename[species],)
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

    def gtcbytaxa(self,taxon):
	accmap = defaultdict(set)
	for acc,gtc in self.allgtc():
	    accmap[acc].add(gtc)
	seen = set()
	for acc,taxid in self.alltaxa():
	    if int(taxid) == int(taxon):
		for gtc in accmap[acc]:
		    if gtc not in seen:
			yield gtc
			seen.add(gtc)

class GlyGenTS(TripleStoreResource):
    endpt = "http://sparql.glygen.org:8880/sparql/query"
    # endpt = "https://sparql.glygen.org/query"
    defns = "http://glygen.org/glycan/"

    def __init__(self,**kw):
        super(GlyGenTS,self).__init__(**kw)
        for k in self.keys():
            self.modify_method(k,partitioner())

    def allglycans(self):
	for row in self.query_glycans():
            yield row['accession']

class GlyGenBetaTS(GlyGenTS):
    endpt = "http://beta-sparql.glygen.org:8880/sparql/query"
    # endpt = "https://beta-sparql.glygen.org/"
    def __init__(self):
        iniFile = os.path.join(os.path.dirname(os.path.realpath(__file__)),"glygen.ini")
        super(GlyGenBetaTS,self).__init__(usecache=False,iniFile=iniFile)

class GlyGen(GlyGenTS):
    pass

class GlyGenBeta(GlyGenBetaTS):
    pass

if __name__ == "__main__":

    cls = sys.argv[1]
    resource = eval(cls+"()")
    query = sys.argv[2]
    method = getattr(resource,query)
    headers = None
    args = sys.argv[3:]
    kwargs = {}
    for i in range(len(args)-1,-1,-1):
	if re.search(r'^[a-z]+=',args[i]):
	    k,v = args[i].split('=',1)
	    kwargs[k] = v
	    del args[i]
    result = method(*args,**kwargs)
    if isinstance(result,basestring) or not hasattr(result,'next'):
	print result
    else:
	for r in result:
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