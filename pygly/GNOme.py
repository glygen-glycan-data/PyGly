
import sys, os, os.path, urllib, rdflib
from collections import defaultdict

class GNOme(object):
    version = "1.0.1"
    referenceowl = "https://raw.githubusercontent.com/glygen-glycan-data/GNOme/V%s/GNOme.owl"%(version,)
    referencefmt = 'xml'

    def __init__(self,resource=None,format=None):
	if not resource:
	    resource = self.referenceowl
	    format = self.referencefmt
	elif not format:
	    format = 'xml'
	self.gnome = rdflib.Graph()
	self.gnome.parse(resource,format=format)
	self.ns = dict()
	self.ns['owl'] = rdflib.Namespace('http://www.w3.org/2002/07/owl#')
	self.ns['gno'] = rdflib.Namespace('http://ontology.glygen.org/gnome/GNO_')
	self.ns['rdf'] = rdflib.Namespace('http://www.w3.org/1999/02/22-rdf-syntax-ns#')
	self.ns['rdfs'] = rdflib.Namespace('http://www.w3.org/2000/01/rdf-schema#')
	self.ns[None] = rdflib.Namespace("")

    def triples(self,subj=None,pred=None,obj=None):
	for (s,p,o) in self.gnome.triples(((self.uri(subj) if subj  else None),(self.uri(pred) if pred  else None),(self.uri(obj) if obj  else None))):
	    yield s,p,o

    def uri(self,id):
	if isinstance(id,rdflib.term.URIRef):
	    return id
	ns,id = id.split(':',1)
	if ns in self.ns:
	    return self.ns[ns][id]
	return self.ns[None][id]

    def accession(self,uri):
	assert uri.startswith(self.ns['gno'])
	return uri[len(self.ns['gno']):]

    def label(self,uri):
	if isinstance(uri,rdflib.term.URIRef):
	    for s,p,o in self.triples(uri,"rdfs:label",None):
	        return unicode(o)
	for ns in self.ns.values():
	    if ns.startswith('http://') and uri.startswith(ns):
		return uri[len(ns):]
	return unicode(uri)

    def nodes(self):
	for s,p,o in self.triples(None,'rdf:type','owl:Class'):
	     acc = self.accession(s)
	     if acc == "00000011":
		continue
             yield acc

    def attributes(self,accession):
	uri = "gno:%s"%(accession,)
        attr = dict()
	for s,p,o in self.triples(uri):
	    plab = self.label(p)
	    if plab != "subClassOf":
	        olab = self.label(o)
	    else:
		olab = self.accession(o)
            attr[plab] = olab
	attr[u'level'] = self.level(accession)
	return attr

    def edges(self):
	for n in self.nodes():
	    for p in self.parents(n):
		yield p,n

    def root(self):
	return "00000001"

    def parents(self,accession):
	uri = "gno:%s"%(accession,)
	for s,p,o in self.triples(uri,'rdfs:subClassOf'):
	    yield self.accession(o)

    def ancestors(self,accession):
	anc = set()
	for p in self.parents(accession):
	    anc.add(p)
	    anc.update(self.ancestors(p))
	return anc

    def children(self,accession):
	uri = "gno:%s"%(accession,)
	for s,p,o in self.triples(None,'rdfs:subClassOf', uri):
	    yield self.accession(s)

    def descendants(self,accession):
	desc = set()
	for c in self.children(accession):
	    desc.add(c)
	    desc.update(self.descendants(c))
	return desc

    def isleaf(self,accession):
	for ch in self.children(accession):
	    return False
	return True

    def isroot(self,accession):
	return accession == self.root()

    def level(self,accession):
	uri = "gno:%s"%(accession,)
	for s,p,o in self.triples(uri,"gno:00000021",None):
	    return " ".join(unicode(self.label(o)).split()[2:])

    def islevel(self,accession,level):
	return self.level(accession) == level

    def ismolecularweight(self,accession):
	return self.islevel(accession,'molecular weight')

    def isbasecomposition(self,accession):
	return self.islevel(accession,'basecomposition')

    def iscomposition(self,accession):
	return self.islevel(accession,'composition')

    def istopology(self,accession):
	return self.islevel(accession,'topology')

    def issaccharide(self,accession):
	return self.islevel(accession,'saccharide')

    def restrict(self,restriction):
	parents = defaultdict(set)
	keep = set()
	for acc in restriction:
	    keep.add(acc)
	    for anc in self.ancestors(acc):
		if self.issaccharide(anc):
		    if anc in restriction:
		        parents[acc].add(anc)
			keep.add(anc)
		else:
		    parents[acc].add(anc)
		    parents[anc].update(self.parents(anc))
		    keep.add(anc)

        toremove = set()
	for n1 in restriction:
	    for n2 in parents[n1]:
		for n3 in parents[n2]:
		    if n3 in parents[n1]:
			toremove.add((n1,n3))
	for n1,n3 in toremove:
	    parents[n1].remove(n3)

	for n in self.nodes():
	    if n not in keep:
		self.gnome.remove((self.uri("gno:"+n),None,None))
	    else:
		self.gnome.remove((self.uri("gno:"+n),self.uri("rdfs:subClassOf"),None))
		for n1 in parents[n]:
		    self.gnome.add((self.uri("gno:"+n),self.uri("rdfs:subClassOf"),self.uri("gno:"+n1)))

    def write(self,handle):
	self.gnome.serialize(handle)

    def dump(self):
        for acc in sorted(g.nodes()):
	    print acc
	    for k,v in sorted(g.attributes(acc).items()):
		print "  %s: %s"%(k,v)

if __name__ == "__main__":
    g = GNOme()
    cmd = sys.argv[1]
    sys.argv.pop(1)
    if cmd == "restrict":

	restriction = set(open(sys.argv[1]).read().split())
	g.restrict(restriction)
	g.write(sys.stdout)
    
    elif cmd == "dump":

	g.dump()

    else:

	print >>sys.stderr, "Bad command: %s"%(cmd,)
	sys.exit(1)
