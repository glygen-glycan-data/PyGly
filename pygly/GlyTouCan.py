#!/bin/env python27

from rdflib import ConjunctiveGraph, Namespace                                                                    
import urllib2, urllib, json
import os.path, sys, traceback, re
import time, random
from collections import defaultdict, Counter
from operator import itemgetter
from hashlib import md5
from StringIO import StringIO
from PIL import Image
from GlycanFormatter import GlycoCTFormat
from Glycan import Glycan


class GlyTouCanCredentialsNotFound(RuntimeError):
    pass

class GlyTouCanRegistrationError(RuntimeError):
    pass

class GlyTouCan(object):
    endpt = 'http://ts.glytoucan.org/sparql'
    api = 'https://api.glytoucan.org/glycan'
    
    def __init__(self,user=None,apikey=None):

	self.g = None
	self.opener = None
	self.user = user
	self.apikey = apikey
	self.delaytime = .2
	self.delaybatch = 1
	self.maxattempts = 3
	self._lastrequesttime = 0
	self._lastrequestcount = 0
	self.alphamap = None

    def setup_sparql(self):
	 self.g = ConjunctiveGraph(store='SPARQLStore')
         self.g.open(self.endpt)

    def setup_api(self,user=None,apikey=None):
	if user == None:
	    user = self.user
	    apikey = self.apikey
	if user == None:
	    user,apikey = self.getcredentials()
        # print user,apikey
	self.password_mgr = urllib2.HTTPPasswordMgrWithDefaultRealm()
	self.password_mgr.add_password(None, self.api, user, apikey)
	self.handler = urllib2.HTTPBasicAuthHandler(self.password_mgr)
	self.opener = urllib2.build_opener(self.handler)
	self.glycoct_format = None

    def _wait(self,delaytime=None):
	if delaytime != None:
	    time.sleep(delaytime)
	    return
	if (self._lastrequestcount % self.delaybatch) == 0 and self._lastrequestcount > 0:
	    time.sleep(self.delaytime)
	self._lastrequesttime = time.time()
	self._lastrequestcount += 1

    def query(self,sparql):
	self._wait()
	if self.g == None:
	    self.setup_sparql()

        attempt = 0
	response = None
	while response == None and attempt < self.maxattempts:
	    try:
		attempt += 1
	        response = self.g.query(sparql)
	    except:
		traceback.print_exc()
	        self._wait(self.delaytime**attempt)

	if response == None:
	    raise IOError("Cannot query SPARQL endpoint")

	return response

    exists_sparql = """
	PREFIX glytoucan: <http://www.glytoucan.org/glyco/owl/glytoucan#>
	
	SELECT ?Saccharide
	WHERE {
   	    ?Saccharide glytoucan:has_primary_id "%(accession)s" .
	}
    """
    def exists(self,accession):
	self._wait()
	if self.g == None:
	    self.setup_sparql()
	response = self.query(self.exists_sparql%dict(accession=accession))
        for row in response.bindings:
	    return True
	return False

    getseq_sparql = """
	PREFIX glycan: <http://purl.jp/bio/12/glyco/glycan#>
	PREFIX glytoucan: <http://www.glytoucan.org/glyco/owl/glytoucan#>
	
	SELECT DISTINCT ?Sequence
	WHERE {
   	    ?Saccharide glytoucan:has_primary_id "%(accession)s" .
   	    ?Saccharide glycan:has_glycosequence ?GlycoSequence .
   	    ?GlycoSequence glycan:has_sequence ?Sequence .
   	    ?GlycoSequence glycan:in_carbohydrate_format glycan:carbohydrate_format_%(format)s .
	}
    """

    def getseq(self,accession,format="wurcs"):
        assert(format in ("wurcs","glycoct","iupac_extended","iupac_condensed"))

	response = self.query(self.getseq_sparql%dict(accession=accession,format=format))
	     
        seqkey = response.vars[0]
	seq = None
        for row in response.bindings:
	    seq = str(row[seqkey].strip())
	    seq = re.sub(r'\n\n+',r'\n',seq)
	    if format == "wurcs" and '~' in seq:
		continue
	    break
	return seq

    getmass_sparql = """
	PREFIX glytoucan: <http://www.glytoucan.org/glyco/owl/glytoucan#>
	
	SELECT ?mass
	WHERE {
   	    ?Saccharide glytoucan:has_primary_id "%(accession)s" .
   	    ?Saccharide glytoucan:has_derivatized_mass ?mass
	}
    """
    def getmass(self,accession):
	response = self.query(self.getmass_sparql%dict(accession=accession))
        masskey = response.vars[0]
        mass = None
        for row in response.bindings:
	    try:
                mass = float(row[masskey].split('/')[-1])
                break
	    except (TypeError,ValueError):
		pass
        return mass

    getmonocount_sparql = """
	PREFIX glycan: <http://purl.jp/bio/12/glyco/glycan#>
	PREFIX glytoucan: <http://www.glytoucan.org/glyco/owl/glytoucan#>
	PREFIX wurcs: <http://www.glycoinfo.org/glyco/owl/wurcs#>
	
	SELECT ?cnt
	WHERE {
   	    ?Saccharide glytoucan:has_primary_id "%(accession)s" .
   	    ?Saccharide glycan:has_glycosequence ?GlycoSequence .
   	    ?GlycoSequence glycan:in_carbohydrate_format glycan:carbohydrate_format_wurcs . 
	    ?GlycoSequence wurcs:RES_count ?cnt
	}
    """
    def getmonocount(self,accession):
	response = self.query(self.getmonocount_sparql%dict(accession=accession))
        key = response.vars[0]
        value = None
        for row in response.bindings:
	    try:
                value = int(row[key])
                break
	    except (TypeError,ValueError):
		pass
        return value

    def getimage(self,accession,notation="cfg",style="extended",avoidcache=False,trials=1):
	assert(notation in ("cfg",) and style in ("compact","normal","extended"))
	self._wait()
n	if trials > 1:
	    avoidcache = True
	imgcnt = Counter()
	hash2img = dict()
	for t in range(trials):
	    rand = ""
	    if avoidcache:
	        rand = "&rand="+("%.8f"%(random.random(),)).split('.',1)[1]
	    try:
	        imgstr = urllib.urlopen("https://glytoucan.org/glycans/%s/image?format=png&notation=%s&style=%s%s"%(accession,notation,style,rand)).read()
	        if len(imgstr) == 0:
	            imgstr = None
	    except IOError:
	        imgstr = None	
	    if imgstr != None:
	        imghash = md5(imgstr).hexdigest().lower()
	        imgcnt[imghash] += 1
		hash2img[imghash] = imgstr
	# print imgcnt
	if len(imgcnt) == 0:
	    return None, None, None
	imgstr = hash2img[max(imgcnt.items(),key=itemgetter(1))[0]]
	pngh = StringIO(imgstr)
	try:
	    pngimg = Image.open(pngh)
        except IOError:
	    return imgstr,None,None
	width,height = pngimg.size
	return imgstr, width, height

    allmotif_sparql = """
	PREFIX glycan: <http://purl.jp/bio/12/glyco/glycan#>
	PREFIX glytoucan: <http://www.glytoucan.org/glyco/owl/glytoucan#>
	PREFIX rdfs: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
	PREFIX rdf: <http://www.w3.org/2000/01/rdf-schema#>
	
	SELECT DISTINCT ?primary_id ?label ?redend
	WHERE {
   	    ?Saccharide glytoucan:has_primary_id ?primary_id . 
	    ?Saccharide rdfs:type glycan:glycan_motif . 
	    ?Saccharide rdf:label ?label .
	    ?Saccharide glytoucan:is_reducing_end ?redend
	}
    """
    def allmotifs(self):
	response = self.query(self.allmotif_sparql)
        key = response.vars[0]
        value = None
        for row in response.bindings:
	    yield tuple(map(str,map(row.get,response.vars)))

    credfile = ".gtccred"
    @staticmethod
    def getcredentials():
	
	# script directory
	dir = os.path.split(sys.argv[0])[0]
        credfile = os.path.join(dir,GlyTouCan.credfile)
	if os.path.exists(credfile):
	    user,apikey = open(credfile).read().split()
	    return user,apikey

	# local directory
        credfile = GlyTouCan.credfile
	if os.path.exists(credfile):
	    user,apikey = open(credfile).read().split()
	    return user,apikey

	# home directory
	dir = os.path.expanduser("~")
        credfile = os.path.join(dir,GlyTouCan.credfile)
	if os.path.exists(credfile):
	    user,apikey = open(credfile).read().split()
	    return user,apikey

	raise GlyTouCanCredentialsNotFound()

    def fixcompwurcs(self,wurcsseq):
	if not self.alphamap:
	    self.alphamap = dict()
	    for i,c in enumerate(range(ord('a'),ord('z')+1)):
    	        self.alphamap[i+1] = chr(c)
    	        self.alphamap[chr(c)] = (i+1)
	    for i,c in enumerate(range(ord('A'),ord('Z')+1)):
                self.alphamap[i+1+26] = chr(c)
                self.alphamap[chr(c)] = (i+1+26)
	prefix,counts,rest = wurcsseq.split('/',2)
        unodes,nodes,edges = counts.split(',')
        nodes = int(nodes)
	assert '0+' in edges
        edges = (nodes-1)
        ambignode = "|".join(map(lambda i: "%s?"%(self.alphamap[i],),range(1,nodes+1)))
        ambigedge = "%s}-{%s"%(ambignode,ambignode)
        ambigedges = [ambigedge]*edges
        return "%s/%s,%d,%d/%s%s"%(prefix,unodes,nodes,edges,rest,"_".join(ambigedges))

    def monotocomp(self,mstr):
	m = re.search(r'-(\d)[abx]_\d-\d',mstr)
        if m:
            mstr = re.sub(r'-\d[abx]_\d-\d','',mstr)
            assert mstr[int(m.group(1))-1] == 'a', mstr
            mstr = list(mstr)
            if m.group(1) == '2':
                mstr[1] = 'U'
            else:
                mstr[0] = 'u'
            mstr = "".join(mstr)
	return mstr

    def monotobasecomp(self,mstr):
	mstr = self.monotocomp(mstr)
	skelplus = mstr.split('_',1)
	skelplus[0] = re.sub(r'[1234]','x',skelplus[0])
	return "_".join(skelplus)

    def wurcscomptrans(self,wurcsseq,monotrans):
	prefix,counts,rest = wurcsseq.split('/',2)
        monos,rest = rest.split(']/',1)
        ids,rest = rest.split('/',1)
        origfreq = defaultdict(int)
        for id in ids.split('-'):
            origfreq[int(id)-1] += 1
        monolist = monos.lstrip('[').split('][')
        # print monolist
        newmonolist = []
	for mstr in monolist:
            newmonolist.append(monotrans(mstr))
	freq = defaultdict(int)
	for i,mstr in enumerate(newmonolist):
	    freq[mstr] += origfreq[i]
	# print freq
	uniq = len(freq)
        total = sum(freq.values())
	newmonolist = sorted(freq,key=lambda k: newmonolist.index(k))
	# print newmonolist
	counts = "%d,%d,0+"%(uniq,total)
	monostr = "".join(map(lambda s: "[%s]"%(s,),newmonolist))
	ids = []
	for i,mstr in enumerate(newmonolist):
	    ids.extend([str(i+1)]*freq[mstr])
	idstr = '-'.join(ids)
	theseq = self.fixcompwurcs("%s/%s/%s/%s/"%(prefix,counts,monostr,idstr))
        # acc,new = self.register(theseq)
	# return self.getseq(acc,'wurcs')
        return theseq

    def makecompwurcs(self,wurcsseq):
	return self.wurcscomptrans(wurcsseq,self.monotocomp)

    def makebasecompwurcs(self,wurcsseq):
	return self.wurcscomptrans(wurcsseq,self.monotobasecomp)

    def haspage(self,accession):
	req = urllib2.Request('https://glytoucan.org/Structures/Glycans/%s'%(accession,))
	self._wait()
	handle = urllib2.urlopen(req)
	# print handle.getcode()
	# print handle.info()
	page = handle.read()
	m = re.search(r'<title>(.*)</title>',page)
	if m and m.group(1).startswith('Accession Number'):
	    return True
	return False

    gettopo_sparql = """
	PREFIX glytoucan: <http://www.glytoucan.org/glyco/owl/glytoucan#>
	PREFIX relation:  <http://www.glycoinfo.org/glyco/owl/relation#>
	
	SELECT ?topo
	WHERE {
   	    ?Saccharide glytoucan:has_primary_id "%(accession)s" .
	    ?Saccharide relation:has_topology ?topo
	}
    """
    def gettopo(self,accession):
        response = self.query(self.gettopo_sparql%dict(accession=accession))
        key = response.vars[0]
        value = None
        for row in response.bindings:
	    value = str(row[key])
            break
	if not value:
	    return None
        return value.rsplit('/',1)[1]

    getcomp_sparql = """
	PREFIX glytoucan: <http://www.glytoucan.org/glyco/owl/glytoucan#>
	PREFIX relation:  <http://www.glycoinfo.org/glyco/owl/relation#>
	
	SELECT ?comp
	WHERE {
   	    ?Saccharide glytoucan:has_primary_id "%(accession)s" .
	    ?Saccharide relation:has_topology ?topo .
	    ?topo relation:has_composition ?comp
	}
    """
    def getcomp(self,accession):
        response = self.query(self.getcomp_sparql%dict(accession=accession))
        key = response.vars[0]
        value = None
        for row in response.bindings:
	    value = str(row[key])
            break
	if not value:
	    return None
        return value.rsplit('/',1)[1]

    def register(self,glycan,user=None,apikey=None):
	if not self.opener:
	    self.setup_api(user=user,apikey=apikey)
	if isinstance(glycan,Glycan):
	    if not self.glycoct_format:
		self.glycoct_format = GlycoCTFormat()
	    sequence = self.glycoct_format.toStr(glycan)
	    # print sequence
	    if not glycan.has_root():
		acc,new = self.register(sequence)
		# print acc,new
		wurcs = self.getseq(acc)
		# print wurcs
		if '0+' in wurcs:
		    sequence = self.fixcompwurcs(wurcs)
		else:
		    return acc,new
	else:
	    sequence = glycan
	sequence = re.sub(r'\n\n+',r'\n',sequence)
	params = json.dumps(dict(sequence=sequence))                                                                  
        # print params
        req = urllib2.Request('https://api.glytoucan.org/glycan/register',params)
    	req.add_header('Content-Type', 'application/json')
        req.add_header('Accept','application/json')
        try:
	    response = None
	    self._wait()
            response = json.loads(self.opener.open(req).read())
	    accession = response['message']
	except (ValueError,IOError), e:
	    # print traceback.format_exc()
	    # print response
	    # force reinitialization of opener...
	    self.opener = None
	    raise GlyTouCanRegistrationError(str(e))
	if response['error'] == "":
	    new = True
	else:
	    new = False
	return accession,new

if __name__ == "__main__":

    import sys

    def items():
        any = False
        for f in sys.argv[1:]:
            any = True
            yield f.strip()
        if not any:
            for l in sys.stdin:
                yield l.strip()

    cmd = sys.argv.pop(1)

    if cmd.lower() == "register":

        gtc = GlyTouCan()
	for f in items():
	    h = open(f); sequence = h.read().strip(); h.close()
	    if not sequence:
		continue
	    try:
		acc,new = gtc.register(sequence)
		print f,acc,("new" if new else "")
	    except GlyTouCanRegistrationError:
		# traceback.print_exc()
		print f,"-","error"
	    sys.stdout.flush()

    elif cmd.lower() in ("wurcs","glycoct"):

        gtc = GlyTouCan()
	for acc in items():
	    print gtc.getseq(acc,cmd)

    elif cmd.lower() in ("image",):

        gtc = GlyTouCan()
	for acc in items():
	    imgstr,width,height = gtc.getimage(acc,trials=5)
	    if imgstr:
	        print acc,width,height
	        wh = open(acc+".png",'w')
		wh.write(imgstr)
		wh.close()

    elif cmd.lower() == "summary":

	gtc = GlyTouCan()
	for acc in items():
	    print "Exists:",gtc.exists(acc)
	    print "WURCS:",bool(gtc.getseq(acc,'wurcs'))
	    print "GlycoCT:",bool(gtc.getseq(acc,'glycoct'))
	    print "Mass:",gtc.getmass(acc)
	    print "Monosaccharide Count:",gtc.getmonocount(acc)
	    print "Composition:",gtc.getcomp(acc)
	    print "Topology:",gtc.gettopo(acc)
	    imgstr,width,height = gtc.getimage(acc,style='extended')
	    print "Extended Image: %s (%dx%d)"%(bool(imgstr),width,height,)
	    imgstr,width,height = gtc.getimage(acc,style='normal')
	    print "Normal Image: %s (%dx%d)"%(bool(imgstr),width,height,)
	    imgstr,width,height = gtc.getimage(acc,style='compact')
	    print "Compact Image: %s (%dx%d)"%(bool(imgstr),width,height,)
	    print "HasPage:",gtc.haspage(acc)

    elif cmd.lower() == "motifs":
	
	gtc = GlyTouCan()
	for acc,label,redend in gtc.allmotifs():
	    print acc,label,redend

    else:
	print >>sys.stderr, "Bad command: %s"%(cmd,)
	sys.exit(1)
	
