
import sys, os, os.path, urllib, rdflib, copy
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

    def get_basecomposition(self,accession):
	bcomps = set()
	for t in self.ancestors(accession):
	    if self.isbasecomposition(t):
		bcomps.add(t)
        exclude = set()
	for t in bcomps:
	    if len(self.descendants(t) & bcomps) > 0:
		exclude.add(t)
	for t in exclude:
	    bcomps.remove(t)
	assert len(bcomps) <= 1
	if len(bcomps) == 0:
	    return None
	return iter(bcomps).next()

    def get_composition(self,accession):
	comps = set()
	for t in self.ancestors(accession):
	    if self.iscomposition(t):
		comps.add(t)
        exclude = set()
	for t in comps:
	    if len(self.descendants(t) & comps) > 0:
		exclude.add(t)
	for t in exclude:
	    comps.remove(t)
	assert len(comps) <= 1
	if len(comps) == 0:
	    return None
	return iter(comps).next()

    def get_topology(self,accession):
	topos = set()
	for t in self.ancestors(accession):
	    if self.istopology(t):
		topos.add(t)
        exclude = set()
	for t in topos:
	    if len(self.descendants(t) & topos) > 0:
		exclude.add(t)
	for t in exclude:
	    topos.remove(t)
	assert len(topos) <= 1
	if len(topos) == 0:
	    return None
	return iter(topos).next()

    def has_basecomposition(self,accession):
	assert self.isbasecomposition(accession)
	return self.descendants(accession)

    def has_composition(self,accession):
	assert self.iscomposition(accession)
	return self.descendants(accession)

    def has_topology(self,accession):
	assert self.istopology(accession)
	return self.descendants(accession)

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

from alignment import GlycanSubsumption, GlycanEqual
from GlyTouCan import GlyTouCan
from Monosaccharide import Anomer
import time

class SubsumptionGraph:
    def __init__(self,*args,**kwargs):
        self.gtc = GlyTouCan(usecache=True)
        self.subsumption = GlycanSubsumption()
        self.geq = GlycanEqual()
        self.verbose = kwargs.get('verbose',0)
        
        masscluster = defaultdict(dict)
        if len(args) > 0:
            for mass in map(float,args):
                rmass = round(mass,2)
                for glyacc in self.gtc.hasmass(rmass):
                    masscluster[rmass][glyacc] = dict(accession=glyacc)
        else:
            for glyacc,mass in self.gtc.allmass():
                rmass = round(mass,2)
                masscluster[rmass][glyacc] = dict(accession=glyacc)
        for rmass,cluster in sorted(masscluster.items()):
            self.compute_component(rmass,cluster)

    def warning(self, msg, level):
        if self.verbose >= level:
            print "# WARNING:%d - %s"%(level,msg)

    def compute_component(self,rmass,cluster):
        start = time.time()

        print "# START %s - %d accessions in molecular weight cluster for %s"%(time.ctime(),len(cluster),rmass)
        sys.stdout.flush()
        
        badparse = 0
        total = len(cluster)
        allgly = dict()
        for acc in cluster:
            gly = self.gtc.getGlycan(acc)
            if not gly:
                badparse += 1
                skels = self.gtc.getUnsupportedSkeletonCodes(acc)
                if len(skels) > 0:
                    for skel in skels:
                        self.warning("unsupported skeleton code: "+skel+" in glycan "+acc,2)
                else:
                    self.warning("unknown problem parsing glycan "+acc,2)
                continue
            cluster[acc]['glycan'] = gly
            cluster[acc]['mass'] = self.gtc.getmass(acc)

        clusteracc = set(map(lambda t: t[0],filter(lambda t: t[1].has_key('glycan'),cluster.items())))

        outedges  = defaultdict(set)
        inedges = defaultdict(set)

        for acc1 in sorted(clusteracc):
            gly1 = cluster[acc1]['glycan']
            for acc2 in sorted(clusteracc):
                gly2 = cluster[acc2]['glycan']
                if acc1 != acc2:
                    if self.subsumption.leq(gly1,gly2):
                        if not self.geq.eq(gly1,gly2) or acc2 < acc1:
                            outedges[acc2].add(acc1)
                            inedges[acc1].add(acc2)
        
        for acc in clusteracc:
            topo = self.gtc.gettopo(acc)
            if topo:
                if topo not in cluster:
                    self.warning("annotated topology %s of %s is not in %s rounded mass cluster"%(topo,acc,rmass),1)
                elif not cluster[topo].get('glycan'):
                    self.warning("annotated topology %s of %s cannot be parsed"%(topo,acc),1)
                elif acc not in outedges[topo] and acc != topo:
                    self.warning("annotated topology %s does not subsume %s"%(topo,acc),1)
            comp = self.gtc.getcomp(acc)
            if comp:
                if comp not in cluster:
                    self.warning("annotated composition %s of %s is not in %s rounded mass cluster"%(comp,acc,rmass),1)
                elif not cluster[comp].get('glycan'):
                    self.warning("annotated composition %s of %s cannot be parsed"%(comp,acc),1)
                elif acc not in outedges[comp] and acc != comp:
                    self.warning("annotated composition %s does not subsume %s"%(comp,acc),1)
            bcomp = self.gtc.getbasecomp(acc)
            if bcomp:
                if bcomp not in cluster:
                    self.warning("annotated base composition %s of %s is not in %s rounded mass cluster"%(bcomp,acc,rmass),1)
                elif not cluster[bcomp].get('glycan'):
                    self.warning("annotated base composition %s of %s cannot be parsed"%(bcomp,acc),1)
                elif acc not in outedges[bcomp] and acc != bcomp:
                    self.warning("annotated base composition %s does not subsume %s"%(bcomp,acc),1)
            try:
                umw = cluster[acc]['glycan'].underivitized_molecular_weight()
            except LookupError:
                umw = None
            if umw == None:
                self.warning("mass could not be computed for %s"%(acc),2)
            elif abs(cluster[acc]['mass'] - umw) > 0.0001:
                self.warning("annotated mass %s for %s is different than computed mass %s"%(cluster[acc]['mass'],acc,umw),1)

        for acc in clusteracc:

            if acc == self.gtc.getbasecomp(acc):
                cluster[acc]['level'] = "BaseComposition"
                cluster[acc]['bcomp'] = acc
            elif acc == self.gtc.getcomp(acc):
                cluster[acc]['level'] = "Composition"
                if self.gtc.getbasecomp(acc):
                    cluster[acc]['bcomp'] = self.gtc.getbasecomp(acc)
                cluster[acc]['comp'] = acc
            elif acc == self.gtc.gettopo(acc):
                cluster[acc]['level'] = "Topology"
                if self.gtc.getbasecomp(acc):
                    cluster[acc]['bcomp'] = self.gtc.getbasecomp(acc)
                if self.gtc.getcomp(acc):
                    cluster[acc]['comp'] = self.gtc.getcomp(acc)
                cluster[acc]['topo'] = acc
            else:
                if self.gtc.gettopo(acc):
                    cluster[acc]['level'] = "Saccharide"
                    cluster[acc]['topo'] = self.gtc.gettopo(acc)
                    if self.gtc.getbasecomp(acc):
                        cluster[acc]['bcomp'] = self.gtc.getbasecomp(acc)
                    if self.gtc.getcomp(acc):
                        cluster[acc]['comp'] = self.gtc.getcomp(acc)

        for acc in clusteracc:
            gly = cluster[acc]['glycan']
            level = cluster[acc].get('level')
            if self.any_anomer(gly):
                if not level:
                    cluster[acc]['level'] = 'Saccharide*'
                elif level != "Saccharide":
                    self.warning("annotation inferred level %s for %s != computed level Saccharide (anomer)"%(level,acc),1)
                continue
            if self.any_parent_pos(gly):
                if not level:
                    cluster[acc]['level'] = 'Saccharide*'
                elif level != "Saccharide":
                    self.warning("annotation inferred level %s for %s != computed level Saccharide (parent_pos)"%(level,acc),1)
                continue
            if self.any_links(gly):
                if not level:
                    cluster[acc]['level'] = 'Topology*'
                elif level != "Topology":
                    self.warning("annotation inferred level %s for %s != computed level Topology"%(level,acc),1)
                continue
            if self.monosaccharide_count(gly) == 1 and self.any_ring(gly):
                if not level:
                    cluster[acc]['level'] = 'Topology*'
                elif level != "Topology":
                    self.warning("annotation inferred level %s for %s != computed level Topology"%(level,acc),1)
                continue
            if self.any_stem(gly):
                if not level:
                    cluster[acc]['level'] = 'Composition*'
                elif level != "Composition":
                    self.warning("annotation inferred level %s for %s != computed level Composition"%(level,acc),1)
                continue
            if not level:
                cluster[acc]['level'] = 'BaseComposition*'
            elif level != "BaseComposition":
                self.warning("annotation inferred level %s for %s != computed level BaseComposition"%(level,acc),1)

        for acc in clusteracc:
            g = cluster[acc]
            if not g.get('topo'):
                if g.get('level') in ("Topology","Topology*"):
                    g['topo'] = acc+"*"
                else:
                    for acc1 in inedges[acc]:
                        if cluster[acc1].get('level') in ("Topology","Topology*"):
                            if not g.get('topo') or g.get('topo').rstrip("*") in inedges[acc1]:
                                g['topo'] = acc1+"*"
            if not g.get('comp'):
                if g.get('level') in ("Composition","Composition*"):
                    g['comp'] = acc+"*"
                else:
                    for acc1 in inedges[acc]:
                        if cluster[acc1].get('level') in ("Composition","Composition*"):
                            if not g.get('comp') or g.get('comp').rstrip("*") in inedges[acc1]:
                                g['comp'] = acc1+"*"
            if not g.get('bcomp'):
                if g.get('level') in ("BaseComposition","BaseComposition*"):
                    g['bcomp'] = acc+"*"
                else:
                    for acc1 in inedges[acc]:
                        if cluster[acc1].get('level') in ("BaseComposition","BaseComposition*"):
                            if not g.get('bcomp') or g.get('bcomp').rstrip("*") in inedges[acc1]:
                                g['bcomp'] = acc1+"*"
        if len(clusteracc) < 1:
            print "# DONE - Elapsed time %.2f sec."%(time.time()-start,)
            sys.stdout.flush()
            return

        print "# NODES - %d/%d glycans in molecular weight cluster for %s"%(len(clusteracc),total,rmass)
        for acc in sorted(clusteracc,key=lambda acc: cluster[acc].get('level')):
            g = cluster[acc]
            print acc,g.get('mass'),g.get('level'),g.get('topo'),g.get('comp'),g.get('bcomp'),
	    gly = g.get('glycan')
            print ("COMP " if (not gly.has_root()) else "") + \
                  ("UNDET " if (gly.undetermined() and gly.has_root()) else "") + \
                  ("FULL" if gly.fully_determined() else "")
        print "# ENDNODES - %d/%d glycans in molecular weight cluster for %s"%(len(clusteracc),total,rmass)
        sys.stdout.flush()

        prunedoutedges = self.prune_edges(outedges)
        prunedinedges  = self.prune_edges(inedges)

        print "# TREE"
        for r in clusteracc:
            if len(prunedinedges[r]) == 0:
                self.print_tree(cluster,prunedoutedges,r)
        print "# ENDTREE"
        
        print "# EDGES"
        for n in clusteracc:
            print "%s:"%(n,),
            print " ".join(prunedoutedges[n])
        print "# ENDEDGES"
        sys.stdout.flush()

        print "# DONE - Elapsed time %.2f sec."%(time.time()-start,)
        sys.stdout.flush()            

    def any_anomer(self,gly):
        for m in gly.all_nodes():
            if m.anomer() in (Anomer.alpha,Anomer.beta):
                return True
        return False

    def any_parent_pos(self,gly):
        for l in gly.all_links():
            if l.parent_pos() != None and l.parent_pos() != set([1]):
                return True
        return False

    def any_ring(self,gly):
        for m in gly.all_nodes():
            if m.ring_start() or m.ring_end():
                return True
        return False

    def any_links(self,gly):
        for l in gly.all_links():
            return True
        return False

    def any_stem(self,gly):
        for m in gly.all_nodes():
            if m.stem():
                return True
        return False

    def monosaccharide_count(self,gly):
	return sum(1 for _ in gly.all_nodes())

    def prune_edges(self,inedges):
        edges = copy.deepcopy(inedges)
        toremove = set()
        for n1 in list(edges):
            for n2 in edges[n1]:
                for n3 in edges[n2]:
                    assert n3 in edges[n1], ", ".join([n1, n2, n3,":"]+edges[n1]+[":"]+edges[n2])
                    if n3 in edges[n1]:
                        toremove.add((n1,n3))
        for n1,n3 in toremove:
            edges[n1].remove(n3)
        return edges

    def print_tree(self,cluster,edges,root,indent=0):
        print "%s%s"%(" "*indent,root),cluster[root]['level']
        for ch in edges[root]:
            self.print_tree(cluster,edges,ch,indent+2)        

if __name__ == "__main__":
    cmd = sys.argv[1]
    sys.argv.pop(1)

    if cmd == "restrict":

        restriction = set(open(sys.argv[1]).read().split())

        g = GNOme()
	g.restrict(restriction)
	g.write(sys.stdout)
    
    elif cmd == "dump":

        g = GNOme()
	g.dump()

    elif cmd == "compute":

        verbose = 0
        while len(sys.argv) > 1 and sys.argv[1] == "-v":
            verbose += 1
            sys.argv.pop(1)

        g = SubsumptionGraph(*sys.argv[1:],verbose=verbose)

    else:

	print >>sys.stderr, "Bad command: %s"%(cmd,)
	sys.exit(1)
