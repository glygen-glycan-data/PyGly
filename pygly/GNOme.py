import copy
import re
import sys
import ssl
import os, os.path, urllib, time
import urllib2
from collections import defaultdict
import datetime

import rdflib
import json

class GNOmeAPI(object):

    # Base class for GNOme and subsumption API supported by both the
    # released OWL file and the so-called Subsumption Graph raw dump
    # file.

    # All nodes
    def nodes(self):
        raise NotImplemented("GNOme API method nodes not implemented")

    # Remove a node
    def delete_node(self, accession):
        raise NotImplemented("GNOme API method delete_node not implemented")

    # Child nodes of a node
    def children(self, accession):
        raise NotImplemented("GNOme API method children not implemented")

    # Parent nodes of a node
    def parents(self, accession):
        raise NotImplemented("GNOme API method parents not implemented")

    # Set the parents of a node
    def set_parents(self, accession, parents):
        raise NotImplemented("GNOme API method set_parents not implemented")

    # Subsumption level of a node
    def level(self, accession):
        raise NotImplemented("GNOme API method level not implemented")

    def get_molecularweight(self, accession):
        raise NotImplemented("GNOme API method get_molecularweight not implemented")

    # Base composition node of a node
    def get_basecomposition(self, accession):
        raise NotImplemented("GNOme API method get_basecomposition not implemented")

    # Composition node of a node
    def get_composition(self, accession):
        raise NotImplemented("GNOme API method get_composition not implemented")

    # Topology node of a node
    def get_topology(self, accession):
        raise NotImplemented("GNOme API method get_topology not implemented")

    ####
    #### Derived functionality, in terms of primitives above...
    ####

    # All edges
    def edges(self):
        for n in self.nodes():
            for c in self.children(n):
                yield n,c

    # Descendants of a node
    def descendants(self, accession):
        desc = set()
        for c in self.children(accession):
            desc.add(c)
            desc.update(self.descendants(c))
        return desc

    # Ancestors of a node
    def ancestors(self, accession):
        anc = set()
        for p in self.parents(accession):
            anc.add(p)
            anc.update(self.ancestors(p))
        return anc

    # Whether a node is a leaf
    def isleaf(self, accession):
        for ch in self.children(accession):
            return False
        return True

    # Whether a node is a root node
    def isroot(self,accession):
        for pt in self.parents(accession):
            return False
        return True

    LEVEL_MOLECULAR_WEIGHT = 'molecular weight'
    LEVEL_BASECOMPOSITION = 'basecomposition'
    LEVEL_COMPOSITION = 'composition'
    LEVEL_TOPOLOGY = 'topology'
    LEVEL_SACCHARIDE = 'saccharide'

    def islevel(self, accession, level):
        return self.level(accession) == level

    def ismolecularweight(self, accession):
        return self.islevel(accession, self.LEVEL_MOLECULAR_WEIGHT)

    def isbasecomposition(self, accession):
        return self.islevel(accession, self.LEVEL_BASECOMPOSITION)

    def iscomposition(self, accession):
        return self.islevel(accession, self.LEVEL_COMPOSITION)

    def istopology(self, accession):
        return self.islevel(accession, self.LEVEL_TOPOLOGY)

    def issaccharide(self, accession):
        return self.islevel(accession, self.LEVEL_SACCHARIDE)

    # All nodes with a base composition node
    def has_molecularweight(self, accession):
        assert self.ismolecularweight(accession)
	if self.get_molecularweight(accession) == accession:
	    yield accession
        for desc in self.descendants(accession):
	    if self.get_molecularweight(desc) == accession:
		yield desc

    # All nodes with a base composition node
    def has_basecomposition(self, accession):
        assert self.isbasecomposition(accession)
	if self.get_basecomposition(accession) == accession:
	    yield accession
        for desc in self.descendants(accession):
	    if self.get_basecomposition(desc) == accession:
		yield desc

    # All nodes with a composition node
    def has_composition(self, accession):
        assert self.iscomposition(accession)
	if self.get_composition(accession) == accession:
	    yield accession
        for desc in self.descendants(accession):
	    if self.get_composition(desc) == accession:
		yield desc

    # All nodes with a topology node
    def has_topology(self, accession):
        assert self.istopology(accession), "Not a topology: "+accession
	if self.get_topology(accession) == accession:
	    yield accession
        for desc in self.descendants(accession):
	    if self.get_topology(desc) == accession:
		yield desc

    def restrict(self, restriction):

        # Update the restriction set to include "landmark" topology, composition, etc nodes
        keep = set(restriction)
        keep.add(self.root())
        for acc in restriction:
            keep.add(self.get_topology(acc))
            keep.add(self.get_composition(acc))
            keep.add(self.get_basecomposition(acc))
            keep.add(self.get_molecularweight(acc))
        if None in keep:
            keep.remove(None)

        # find all ancestors of each kept node
        parents = defaultdict(set)
        for acc in keep:
            for anc in self.ancestors(acc):
		if anc in keep:
                    parents[acc].add(anc)

        # then eliminate shortcuts
        toremove = set()
        for n1 in keep:
            for n2 in parents[n1]:
                for n3 in parents[n2]:
                    if n3 in parents[n1]:
                        toremove.add((n1, n3))
        for n1, n3 in toremove:
            parents[n1].remove(n3)

        for n in self.nodes():
            if n not in keep:
                self.delete_node(n)
            else:
                self.set_parents(n,parents[n])

class GNOme(GNOmeAPI):
    referenceowl = "http://purl.obolibrary.org/obo/gno.owl"
    # referenceowl = "http://github.com/glygen-glycan-data/GNOme/releases/latest/download/GNOme.owl"
    referencefmt = 'xml'

    def __init__(self, resource=None, format=None):
        if not resource:
            resource = self.referenceowl
            format = self.referencefmt
        elif not format:
            format = 'xml'
        self.gnome = rdflib.Graph()
        self.gnome.parse(resource, format=format)
        self.ns = dict()
        self.ns['owl'] = rdflib.Namespace('http://www.w3.org/2002/07/owl#')
        self.ns['gno'] = rdflib.Namespace('http://purl.obolibrary.org/obo/GNO_')
        self.ns['rdf'] = rdflib.Namespace('http://www.w3.org/1999/02/22-rdf-syntax-ns#')
        self.ns['rdfs'] = rdflib.Namespace('http://www.w3.org/2000/01/rdf-schema#')
        self.ns[None] = rdflib.Namespace("")
	versionurl = None
        for s,p,o in self.triples("owl:Ontology","owl:versionIRI",None):
            versionurl = str(o)
            break
	if versionurl:
            self.version = versionurl.split('/')[-2][1:]

    def triples(self, subj=None, pred=None, obj=None):
        for (s, p, o) in self.gnome.triples(((self.uri(subj) if subj else None),
                                             (self.uri(pred) if pred else None),
                                             (self.uri(obj) if obj else None))):
            yield s, p, o

    def uri(self, id):
        if isinstance(id, rdflib.term.URIRef):
            return id
        ns, id = id.split(':', 1)
        if ns in self.ns:
            return self.ns[ns][id]
        return self.ns[None][id]

    def accession(self, uri):
        assert uri.startswith(self.ns['gno'])
        return uri[len(self.ns['gno']):]

    def label(self, uri):
        if isinstance(uri, rdflib.term.URIRef):
            for s, p, o in self.triples(uri, "rdfs:label", None):
                return unicode(o)
        for ns in self.ns.values():
            if ns.startswith('http://') and uri.startswith(ns):
                return uri[len(ns):]
        return unicode(uri)

    def nodes(self):
        for s, p, o in self.triples(None, 'rdf:type', 'owl:Class'):
            acc = self.accession(s)
            if acc == "00000011":
                continue
            yield acc

    def attributes(self, accession):
        uri = "gno:%s" % (accession,)
        attr = dict()

        for s, p, o in self.triples(uri):
            plab = self.label(p)
            if plab != "subClassOf":
                olab = self.label(o)
                attr[plab] = olab
            else:
                olab = self.accession(o)
                if plab not in attr:
		                attr[plab] = []
                attr[plab].append(olab)
        attr[u'level'] = self.level(accession)
        return attr

    def root(self):
        return "00000001"

    def parents(self, accession):
        uri = "gno:%s" % (accession,)
        for s, p, o in self.triples(uri, 'rdfs:subClassOf', None):
            yield self.accession(o)

    def children(self, accession):
        uri = "gno:%s" % (accession,)
        for s, p, o in self.triples(None, 'rdfs:subClassOf', uri):
            yield self.accession(s)

    def level(self, accession):
        uri = "gno:%s" % (accession,)
        for s, p, o in self.triples(uri, "gno:00000021", None):
            return self.label(o)

    def get_molecularweight(self, accession):
        # There is a single molecular weight node at or above any node
        # (that is not the root), so it can be determined by the hierarchy
        if self.ismolecularweight(accession):
            return accession
        for anc in self.ancestors(accession):
            if self.ismolecularweight(anc):
                return anc

    def get_basecomposition(self, accession):
        if self.isbasecomposition(accession):
            return accession
        uri = "gno:%s" % (accession,)
        for s, p, o in self.triples(uri, "gno:00000033", None):
            return self.label(o)

    def get_composition(self, accession):
        if self.iscomposition(accession):
            return accession
        uri = "gno:%s" % (accession,)
        for s, p, o in self.triples(uri, "gno:00000034", None):
            return self.label(o)

    def get_topology(self, accession):
        if self.istopology(accession):
            return accession
        uri = "gno:%s" % (accession,)
        for s, p, o in self.triples(uri, "gno:00000035", None):
            return self.label(o)

    def delete_node(self, accession):
        self.gnome.remove((self.uri("gno:" + accession), None, None))

    def set_parents(self, accession, parents):
        self.gnome.remove((self.uri("gno:" + accession), self.uri("rdfs:subClassOf"), None))
        for parent in parents:
            self.gnome.add((self.uri("gno:" + accession), self.uri("rdfs:subClassOf"), self.uri("gno:" + parent)))

    def write(self, handle):
        writer = OWLWriter()
        writer.write(handle, self.gnome)

    def dump(self):
        for acc in sorted(g.nodes()):
            print acc
            for k, v in sorted(g.attributes(acc).items()):
                print "  %s: %s" % (k, v)

    def get_cb_button_str(self, acc):
        for s, p, o in self.triples(None, "gno:00000101", None):
            if acc == self.accession(s):
                return o

    # ics2dp: IUPAC composition string to dictionary pattern(regex)
    ics2dp = re.compile(r"(\D{1,8})(\d{1,3})")
    def iupac_composition_str_to_dict(self, s):
        res = {}
        for p in self.ics2dp.findall(s):
            res[p[0]] = int(p[1])
        return res

    # cbbutton -> composition browser button
    def get_cbbutton(self, acc):
        return self.iupac_composition_str_to_dict(self.get_cb_button_str(acc))

    def all_cbbutton(self):
        res = {}
        for s, p, o in self.triples(None, "gno:00000101", None):
            acc = self.accession(s)
            res[acc] = self.iupac_composition_str_to_dict(o)
        return res

    def toViewerData(self, output_file_path):
        res = {}
        cbbutton = self.all_cbbutton()
        # byonic = self.all_Byonic()
        syms = self.all_synonym()

        for n in self.nodes():
            t = ""
            if self.isbasecomposition(n):
                t = "basecomposition"
            elif self.iscomposition(n):
                t = "composition"
            elif self.istopology(n):
                t = "topology"
            elif self.issaccharide(n):
                t = "saccharide"
            else:
                continue

            if n not in cbbutton:
                continue

            children = list(self.children(n))
            res[n] = {
                "level": t, "children": children, "count": cbbutton[n]
            }
            if len(children) == 0:
                del res[n]['children']

            if n in syms:
                res[n]['syms'] = syms[n]

        json.dump(res, open(output_file_path, 'w'))
        return

    def get_Byonic(self, acc):
        for s, p, o in self.triples('gno:'+acc, "gno:00000202", None):
            yield o

    def all_Byonic(self):
        res = {}
        for s, p, o in self.triples(None, "gno:00000202", None):
            acc = self.accession(s)
            res[acc] = o
        return res

    def get_synonym(self, acc):
        for s, p, o in self.triples('gno:'+acc, rdflib.URIRef("http://www.geneontology.org/formats/oboInOwl#hasExactSynonym"), None):
            yield o

    def all_synonym(self):
        res = defaultdict(list)
        for s, p, o in self.triples(None, rdflib.URIRef("http://www.geneontology.org/formats/oboInOwl#hasExactSynonym"), None):
            acc = self.accession(s)
            res[acc].append(o)
        res = dict(res)
        return res




from alignment import GlycanSubsumption, GlycanEqual, GlycanEqualWithWURCSCheck
from Monosaccharide import Anomer
from GlycanResource import GlyTouCan
from manipulation import Topology, Composition, BaseComposition


class SubsumptionGraph(GNOmeAPI):
    def __init__(self, *args, **kwargs):
        pass

    def compute(self, *args, **kwargs):
        self.gtc = GlyTouCan(usecache=False)
        self.subsumption = GlycanSubsumption()
        self.geq = GlycanEqual()
        self.geqwwc = GlycanEqualWithWURCSCheck()
        self.topology = Topology()
        self.composition = Composition()
        self.basecomposition = BaseComposition()
        self.verbose = kwargs.get('verbose', 0)
        self.score = IncompleteScore()

	# invalid = set(self.gtc.allinvalid())
	# invalid = set()
	replace = dict(self.gtc.replace())

        masscluster = defaultdict(dict)
	clustermap = defaultdict(set)
        if len(args) > 0:
	    argmass = []
	    for a in args:
		if '-' in a:
		    low,high = [ s.strip() for s in split('-') ]
		    if not low:
			low = 0.0
		    if not high:
			high = 1e+20
		    low = str(round(float(low),2))
		    high = str(round(float(high),2))
		else:
		    low = str(round(float(a),2))
		    high = low
		argmass.append((low,high))
	    for glyacc in self.gtc.allaccessions():
		mass = self.gtc.getmass(glyacc)
		if not mass:
		    mass = self.gtc.umw(glyacc,fetch='wurcs')
	        if not mass:
		    continue
		rmass = str(round(mass, 2))
		if glyacc in clustermap:
		    for rmi in clustermap[glyacc]:
		        if glyacc in masscluster[rmi]:
                            del masscluster[rmi][glyacc]
		    clustermap[glyacc].add(rmass)
		    continue
		for low,high in argmass:
		    if float(low) <= float(rmass) <= float(high):
                        masscluster[rmass][glyacc] = dict(accession=glyacc)
			clustermap[glyacc].add(rmass)
			break
        else:
	    for glyacc in self.gtc.allaccessions():
		mass = self.gtc.getmass(glyacc)
		if not mass:
		    mass = self.gtc.umw(glyacc,fetch='wurcs')
		if not mass:
		    continue
                rmass = str(round(mass, 2))
		if glyacc in clustermap:
		    for rmi in clustermap[glyacc]:
                        if glyacc in masscluster[rmi]:
                            del masscluster[rmi][glyacc]                                                                    
                    clustermap[glyacc].add(rmass)
                    continue
                masscluster[rmass][glyacc] = dict(accession=glyacc)
		clustermap[glyacc].add(rmass)

	for rmass,cluster in masscluster.items():
	    for acc in list(cluster):
		topo = self.gtc.gettopo(acc)
		if topo in replace:
		    topo = replace[topo]
	        if topo and topo not in cluster:
		    masscluster[rmass][topo] = dict(accession=topo)
		    for rmi in clustermap.get(topo,[]):
			del masscluster[rmi][topo]
		comp = self.gtc.getcomp(acc)
		if comp in replace:
		    comp = replace[comp]
	        if comp and comp not in cluster:
		    masscluster[rmass][comp] = dict(accession=comp)
		    for rmi in clustermap.get(comp,[]):
		        del masscluster[rmi][comp]
		bcomp = self.gtc.getbasecomp(acc)
		if bcomp in replace:
		    bcomp = replace[bcomp]
	        if bcomp and bcomp not in cluster:
		    masscluster[rmass][bcomp] = dict(accession=bcomp)
		    for rmi in clustermap.get(bcomp,[]):
		        del masscluster[rmi][bcomp]

        for rmass, cluster in sorted(masscluster.items(),key=lambda t: float(t[0])):
            self.compute_component(rmass, cluster)

    def warning(self, msg, level):
        if self.verbose >= level:
            print "# WARNING:%d - %s" % (level, msg)
	    sys.stdout.flush()

    def compute_component(self, rmass, cluster):
        start = time.time()

        replace = self.gtc.replace()

        print "# START %s - %d accessions in molecular weight cluster for %s" % (time.ctime(), len(cluster), rmass)
        sys.stdout.flush()

        badparse = 0
        total = len(cluster)
        allgly = dict()
        for acc in sorted(cluster):
            gly = self.gtc.getGlycan(acc,format='wurcs')
            # gly = self.gtc.getGlycan(acc,format='glycoct')
            if not gly:
                badparse += 1
                skels, substs, invalid, other = self.gtc.getUnsupportedCodes(acc)
                if len(skels) > 0:
                    for skel in skels:
                        self.warning("unsupported skeleton code: " + skel + " in glycan " + acc, 2)
		if len(substs) > 0:
		    for subst in substs:
                        self.warning("unsupported substituent: " + subst + " in glycan " + acc, 2)
		if len(invalid) > 0:
		    for inv in invalid:
                        self.warning("unsupported monosaccharide: " + inv + " in glycan " + acc, 2)
		if len(other) > 0:
		    for oth in other:
                        self.warning("other glycan error: " + oth + " in glycan " + acc, 2)
		if max(map(len,[skels,substs,invalid,other])) == 0:
                    self.warning("unknown problem parsing glycan " + acc, 2)
                continue
            cluster[acc]['glycan'] = gly
	    mass = self.gtc.getmass(acc)
	    mymass = self.gtc.umw(acc)
	    if mass and abs(float(rmass)-mass) <= 0.01:
                cluster[acc]['mass'] = mass
	    elif mymass and abs(float(rmass)-mymass) <= 0.01:
                cluster[acc]['mass'] = mymass
	    else:
		cluster[acc]['mass'] = float(rmass)

        clusteracc = set(map(lambda t: t[0], filter(lambda t: t[1].has_key('glycan'), cluster.items())))

        outedges = defaultdict(set)
        inedges = defaultdict(set)

	self.warning("Computation of subsumption relationships started",5)
        for acc1 in sorted(clusteracc):
            gly1 = cluster[acc1]['glycan']
            for acc2 in sorted(clusteracc):
                gly2 = cluster[acc2]['glycan']
                if acc1 != acc2:
		    self.warning("%s <?= %s"%(acc1,acc2),5)
                    if self.subsumption.leq(gly1, gly2):
			iseq = self.geq.eq(gly1, gly2)
                        if not iseq or acc2 < acc1:
                            outedges[acc2].add(acc1)
                            inedges[acc1].add(acc2)
                        if iseq and self.geqwwc.eq(gly1, gly2) and acc1 < acc2:
                            self.warning("Potential WURCS canonicalization issue: %s == %s"%(acc1,acc2),1)

	self.warning("Computation of subsumption relationships done",5)

        # Check GlyTouCan topology, composition, basecomposition, and
        # molecular weight annotations with respect to computed
        # subsumption relationships and computed molecular weight

        for acc in sorted(clusteracc):

            topo = self.gtc.gettopo(acc)
	    if topo in replace:
		topo = replace[topo]
            if topo:
                if topo not in cluster:
                    self.warning("annotated topology %s of %s is not in %s rounded mass cluster" % (topo, acc, rmass), 1)
                elif not cluster[topo].get('glycan'):
                    self.warning("annotated topology %s of %s cannot be parsed" % (topo, acc), 1)
                elif acc not in outedges[topo] and acc != topo:
                    self.warning("annotated topology %s does not subsume %s" % (topo, acc), 1)

            comp = self.gtc.getcomp(acc)
	    if comp in replace:
		comp = replace[comp]
            if comp:
                if comp not in cluster:
                    self.warning("annotated composition %s of %s is not in %s rounded mass cluster" % (comp, acc, rmass), 1)
                elif not cluster[comp].get('glycan'):
                    self.warning("annotated composition %s of %s cannot be parsed" % (comp, acc), 1)
                elif acc not in outedges[comp] and acc != comp:
                    self.warning("annotated composition %s does not subsume %s" % (comp, acc), 1)

            bcomp = self.gtc.getbasecomp(acc)
	    if bcomp in replace:
		bcomp = replace[bcomp]
            if bcomp:
                if bcomp not in cluster:
                    self.warning("annotated base composition %s of %s is not in %s rounded mass cluster" % (bcomp, acc, rmass),1)
                elif not cluster[bcomp].get('glycan'):
                    self.warning("annotated base composition %s of %s cannot be parsed" % (bcomp, acc), 1)
                elif acc not in outedges[bcomp] and acc != bcomp:
                    self.warning("annotated base composition %s does not subsume %s" % (bcomp, acc), 1)

            try:
                umw = cluster[acc]['glycan'].underivitized_molecular_weight()
            except LookupError:
                umw = None
            if umw == None:
                self.warning("mass could not be computed for %s" % (acc), 2)
            elif cluster[acc]['mass'] == None or abs(cluster[acc]['mass'] - umw) > 0.0001:
                self.warning("annotated mass %s for %s is different than computed mass %s" % (cluster[acc]['mass'], acc, umw), 1)

        # Infer subsumption level and and store topology, composition,
        # basecomposition based on GlyTouCan annotations.

        for acc in sorted(clusteracc):

	    bcomp = self.gtc.getbasecomp(acc)
	    if bcomp in replace:
		bcomp = replace[bcomp]
	    comp = self.gtc.getcomp(acc)
	    if comp in replace:
		comp = replace[comp]
	    topo = self.gtc.gettopo(acc)
	    if topo in replace:
		topo = replace[topo]

            if acc == bcomp:
                cluster[acc]['level'] = "BaseComposition"
                cluster[acc]['bcomp'] = acc

            elif acc == comp:
                cluster[acc]['level'] = "Composition"
                if bcomp:
                    cluster[acc]['bcomp'] = bcomp
                cluster[acc]['comp'] = acc

            elif acc == topo:
                cluster[acc]['level'] = "Topology"
                if bcomp:
                    cluster[acc]['bcomp'] = bcomp
                if comp:
                    cluster[acc]['comp'] = comp
                cluster[acc]['topo'] = acc

            else:
                if topo:
                    cluster[acc]['level'] = "Saccharide"
                    cluster[acc]['topo'] = topo
                    if bcomp:
                        cluster[acc]['bcomp'] = bcomp
                    if comp:
                        cluster[acc]['comp'] = comp

        # Augment GlyTouCan level annotations with levels inferred
        # from glycan structure characteristics. Check to see these
        # are consistent with levels inferred from GlyTouCan
        # annotations. Put a * on those levels we set. Currently,
        # GlyTouCan annotations take precedence.

        for acc in sorted(clusteracc):
            gly = cluster[acc]['glycan']
            level = cluster[acc].get('level')

            if self.any_anomer(gly):
                if not level:
                    cluster[acc]['level'] = 'Saccharide*'
                elif level != "Saccharide":
                    self.warning("annotation inferred level %s for %s != computed level Saccharide (anomer)" % (level, acc), 1)
		    # cluster[acc]['level'] = 'Saccharide*'
                continue

            if self.any_parent_pos(gly):
                if not level:
                    cluster[acc]['level'] = 'Saccharide*'
                elif level != "Saccharide":
                    self.warning("annotation inferred level %s for %s != computed level Saccharide (parent_pos)" % (level, acc),
                        1)
                    # cluster[acc]['level'] = 'Saccharide*'
                continue

            if self.any_links(gly):
                if not level:
                    cluster[acc]['level'] = 'Topology*'
                elif level != "Topology":
                    self.warning("annotation inferred level %s for %s != computed level Topology" % (level, acc), 1)
                    # cluster[acc]['level'] = 'Topology*'
                continue

            if self.mono_count(gly) == 1 and self.any_ring(gly):
                if not level:
                    cluster[acc]['level'] = 'Topology*'
                elif level != "Topology":
                    self.warning("annotation inferred level %s for %s != computed level Topology" % (level, acc), 1)
                    # cluster[acc]['level'] = 'Topology*'
                continue

	    if self.any_ring(gly):
                self.warning("%s has no linkages but does have ring values" % (acc,), 1)

            if self.any_stem(gly):
                if not level:
                    cluster[acc]['level'] = 'Composition*'
                elif level != "Composition":
                    self.warning("annotation inferred level %s for %s != computed level Composition" % (level, acc), 1)
                    # cluster[acc]['level'] = 'Composition*'
                continue

            if not level:
                cluster[acc]['level'] = 'BaseComposition*'
            elif level != "BaseComposition":
                self.warning("annotation inferred level %s for %s != computed level BaseComposition" % (level, acc), 1)
                # cluster[acc]['level'] = 'BaseComposition*'

        # Add topology, composition, basecomposition annotations for
        # those missing from GlyTouCan. Check GlyTouCan annotations to
        # make sure they are consistent with ours. Put a * on those we
        # set. Currently, GlyTouCan annotations take precedence.

	self.warning("Checking topo, comp, bcomp relationships...",5)
        for acc in sorted(clusteracc):
	    self.warning("Checking topo, comp, bcomp relationships for %s"%(acc,),5)
            g = cluster[acc]

            gly = g['glycan']

            topo = set()
            comp = set()
            bcomp = set()
            for acc1 in (list(inedges[acc]) + [acc]):
                level1 = cluster[acc1].get('level').strip("*")
                if level1 == "Saccharide":
                    continue
                gly1 = cluster[acc1]['glycan']

                if False and level1 == "Topology" and acc == "G65022XW":
                    print self.geq.eq(gly1,self.topology(gly))
                    print "topology(%s)"%(acc,)
                    print self.topology(gly).glycoct()
                    print acc1
                    print gly1.glycoct()

                if False and level1 == "Composition" and acc == "G76453ML":
                    print self.geq.eq(gly1,self.composition(gly))
                    print "composition(%s)"%(acc,)
                    print self.composition(gly).glycoct()
                    print acc1
                    print gly1.glycoct()
		
	        self.warning("Checking topo, comp, bcomp relationships for %s: %s (%s)"%(acc,acc1,level1),5)
                if level1 == "BaseComposition" and self.geq.eq(gly1,self.basecomposition(gly)) and acc1 not in replace:
                    bcomp.add(acc1)
                elif level1 == "Composition" and self.geq.eq(gly1,self.composition(gly)) and acc1 not in replace:
                    comp.add(acc1)
                elif level1 == "Topology" and self.geq.eq(gly1,self.topology(gly)) and acc1 not in replace:
                    topo.add(acc1)

            if len(topo) > 1:
                self.warning("multiple topologies %s for %s"%(", ".join(topo), acc), 1)
                if g.get('topo') in topo:
		    # Take GTC one if present
                    topo = g.get('topo')
                else:
		    # Take the one subsumed by all the others...
		    topo1=None
		    for acc in topo:
			if (topo-inedges[acc]) == set([acc]):
			    topo1 = acc
			    break
		    topo = topo1
            elif len(topo) == 0:
                topo = None
            else:
                topo = topo.pop()

            if len(comp) > 1:
                self.warning("multiple compositions %s for %s"%(", ".join(comp), acc), 1)
                if g.get('comp') in comp:
		    # Take GTC one if present
                    comp = g.get('comp')
                else:
		    comp1=None
		    for acc in comp:
			if (comp-inedges[acc]) == set([acc]):
			    comp1 = acc
			    break
                    comp = comp1
            elif len(comp) == 0:
                comp = None
            else:
                comp = comp.pop()

            if len(bcomp) > 1:
                self.warning("multiple base compositions %s for %s"%(", ".join(bcomp), acc), 1)
                if g.get('bcomp') in bcomp:
		    # Take GTC one if present
                    bcomp = g.get('bcomp')
                else:
		    bcomp1=None
		    for acc in bcomp:
			if (bcomp-inedges[acc]) == set([acc]):
			    bcomp1 = acc
			    break
                    bcomp = bcomp1
            elif len(bcomp) == 0:
                bcomp = None
            else:
                bcomp = bcomp.pop()

            if g.get('topo') and g.get('topo') != topo and g.get('topo') in cluster:
                self.warning("annotated topology %s for %s != computed topology %s"%(g.get('topo'),acc,topo),1)
            if acc == topo and cluster[acc].get('level').strip("*") != "Topology":
                self.warning("bad level %s for %s with topology %s"%(cluster[acc].get('level'),acc,topo),1)
            if not g.get('topo') and topo:
                g['topo'] = topo + "*"

            if g.get('comp') and g.get('comp') != comp:
                self.warning("annotated composition %s for %s != computed composition %s"%(g.get('comp'),acc,comp),1)
	        if g.get('comp') in cluster and 'glycan' in cluster[g.get('comp')]:
		    glycomp = cluster[g.get('comp')]['glycan']
		    floating = set()
		    for r in glycomp.all_nodes(undet_subst=True):
		        if not r.is_monosaccharide():
			    if r.name() not in Composition.floating_substs:
			        floating.add(str(r).split(":")[1])
		    assert len(floating) == 0, "Unaccounted for floating substs: "+", ".join(floating)
            if acc == comp and cluster[acc].get('level').strip("*") != "Composition":
                self.warning("bad level %s for %s with composition %s"%(cluster[acc].get('level'),acc,comp),1)
            if not g.get('comp') and comp:
                g['comp'] = comp + "*"

            if g.get('bcomp') and g.get('bcomp') != bcomp and g.get('bcomp') in cluster:
                self.warning("annotated base composition %s for %s != computed base composition %s"%(g.get('bcomp'),acc,bcomp),1)
            if acc == bcomp and cluster[acc].get('level').strip("*") != "BaseComposition":
                self.warning("bad level %s for %s with topology %s"%(cluster[acc].get('level'),acc,bcomp),1)
            if not g.get('bcomp') and bcomp:
                g['bcomp'] = bcomp + "*"

	self.warning("Checking topo, comp, bcomp relationships... done.",5)


        for acc, content in cluster.items():
            try:
                content['missing'] = self.score.score(content["glycan"])
            except:
                content['missing'] = None
                self.warning("Unable to compute missing rank for %s" % (acc), 1)

        if len(clusteracc) < 1:
            print "# DONE - Elapsed time %.2f sec." % (time.time() - start,)
            sys.stdout.flush()
            return

        print "# NODES - %d/%d glycans in molecular weight cluster for %s" % (len(clusteracc), total, rmass)
        for acc in sorted(clusteracc, key=lambda acc: (cluster[acc].get('level').rstrip('*'),acc) ):
            g = cluster[acc]
            print acc,
	    print g.get('mass'),
            print g.get('level'), g.get('topo'), g.get('comp'), g.get('bcomp'),
            gly = g.get('glycan')
	    extras = []
	    if not gly.has_root():
		extras.append("COMP")
	    if gly.undetermined() and gly.has_root():
	        extras.append("UNDET")
	    if gly.fully_determined():
		extras.append("FULL")
	    for m,c in sorted(gly.iupac_composition(floating_substituents=False,aggregate_basecomposition=True).items()):
                if m != "Count" and c > 0:
		    extras.append("%s:%d"%(m,c))
	    print " ".join(extras)
        print "# ENDNODES - %d/%d glycans in molecular weight cluster for %s" % (len(clusteracc), total, rmass)
        sys.stdout.flush()

        prunedoutedges = self.prune_edges(outedges)
        prunedinedges = self.prune_edges(inedges)

        print "# TREE"
        for r in sorted(clusteracc):
            if len(prunedinedges[r]) == 0:
                self.print_tree(cluster, prunedoutedges, r)
        print "# ENDTREE"

        print "# EDGES"
        for n in sorted(clusteracc):
            print "%s:" % (n,),
            print " ".join(sorted(prunedoutedges[n]))
        print "# ENDEDGES"
        sys.stdout.flush()

        print "# DONE - Elapsed time %.2f sec." % (time.time() - start,)
        sys.stdout.flush()

    def any_anomer(self, gly):
        for m in gly.all_nodes():
            if m.anomer() in (Anomer.alpha, Anomer.beta):
                return True
        return False

    def any_parent_pos(self, gly):
        for l in gly.all_links():
            if l.parent().is_monosaccharide():
                pass
            else:
                continue
                l = l.parent().parent_links()[0]
            if l.parent_pos() != None and l.parent_pos() != set([l.parent().ring_start()]):
                return True
        return False

    def any_ring(self, gly):
        for m in gly.all_nodes():
            if m.ring_start() != None or m.ring_end() != None:
                return True
        return False

    def any_links(self, gly):
        for l in gly.all_links():
            return True
        return False

    def any_stem(self, gly):
        for m in gly.all_nodes():
            if m.stem():
                return True
        return False

    def mono_count(self, gly):
        return sum(1 for _ in gly.all_nodes())

    def prune_edges(self, inedges):
        edges = copy.deepcopy(inedges)
        toremove = set()
        for n1 in list(edges):
            for n2 in edges[n1]:
                for n3 in edges[n2]:
                    assert n3 in edges[n1], ", ".join([n1, n2, n3, ":"] + list(edges[n1]) + [":"] + list(edges[n2]))
                    if n3 in edges[n1]:
                        toremove.add((n1, n3))
        for n1, n3 in toremove:
            edges[n1].remove(n3)
        return edges

    def print_tree(self, cluster, edges, root, indent=0):
        print "%s%s" % (" " * indent, root), cluster[root]['level'], cluster[root]['missing']
        for ch in sorted(edges[root]):
            self.print_tree(cluster, edges, ch, indent + 2)


    # Dump file parsing part starts here
    raw_data = {}
    dumpfilepath = ""
    warnings = defaultdict(set)
    warningsbytype = None
    monosaccharide_count = {}
    monosaccharide_count_pattern = re.compile(r"\w{1,8}:\d{1,3}")

    def loaddata(self, dumpfilepath):
        self.readfile(dumpfilepath)
        self.dumpfilepath = dumpfilepath
        for index in range(len(self.errorcategory)):
            self.errorcategory[index][1] = re.compile(self.errorcategory[index][1])

    def readfile(self, dumpfilepath):
        f = open(dumpfilepath)
        raw_data = {}
        for l in f:
            l = l.strip()

            if l.startswith("#"):
                if l.startswith("# WARNING"):
                    t1, msg = l.split(" - ")
                    temp1 = re.compile(r"\d").findall(t1)
                    if len(temp1) != 1:
                        raise RuntimeError
                    level = temp1[0]
                    self.warnings[level].add(msg)
                    continue

                lineInfo = "node"  # node type none
                if l.startswith("# NODES"):
                    lineInfo = "node"
                    mass = float(re.compile("\d{2,5}\.\d{1,2}").findall(l)[0])
                    mass = "%.2f" % mass
                    content = {"nodes": {}, "edges": {}, "missing": {}}
                    raw_data[mass] = content
                elif l.startswith("# EDGES"):
                    lineInfo = "edge"
                elif l.startswith("# TREE"):
                    lineInfo = "tree"
                else:
                    lineInfo = "none"
            else:

                if lineInfo == "none":
                    continue
                elif lineInfo == "node":
                    nodeacc = l.split()[0]
		    nodemw = "%.2f"%(float(l.split()[1]),)
                    nodetype = l.split()[2].rstrip("*")
                    topoacc = l.split()[3].rstrip("*")
		    if topoacc == "None":
			topoacc = None
                    compacc = l.split()[4].rstrip("*")
		    if compacc == "None":
			compacc = None
                    bcompacc = l.split()[5].rstrip("*")
		    if bcompacc == "None":
			bcompacc = None
                    content["nodes"][nodeacc] = (nodemw,nodetype,topoacc,compacc,bcompacc)

                    mono_count = {}
                    for cell in l.split():
                        temp2 = self.monosaccharide_count_pattern.findall(cell)
                        if len(temp2) == 1:
                            temp2 = temp2[0].split(":")
                            iupac_comp_mono, count = temp2[0], int(temp2[1])
                            mono_count[iupac_comp_mono] = count
                        if mono_count:
                            temp3 = ""
                            for iupac_mono_str, mono_count_eatch in sorted(mono_count.items()):
                                temp3 += iupac_mono_str + str(mono_count_eatch)
                            self.monosaccharide_count[nodeacc] = mono_count

                elif lineInfo == "edge":
                    to = re.compile("G\d{5}\w{2}").findall(l)
                    fromx = to.pop(0)
                    if to:
                        content["edges"][fromx] = to

                elif lineInfo == "tree":
                    gtcacc, g_type, asterik, missingrank = list(re.compile("(G\d{5}\w{2}) (BaseComposition|Composition|Topology|Saccharide)(\*)? (\d{1,5})").findall(l))[0]
                    content["missing"][gtcacc] = missingrank

                else:
                    raise RuntimeError
        self.raw_data = raw_data

        self.allnodestype = {}
        self.alledges = {}
        self.allmissingrank = {}
        for component in raw_data.values():
            self.allnodestype.update(copy.deepcopy(component["nodes"]))
            self.alledges.update(copy.deepcopy(component["edges"]))
            self.allmissingrank.update(copy.deepcopy(component["missing"]))

        allmass = list()
        for m in self.raw_data.keys():
            self.allnodestype[m] =  (m, "molecular weight", None, None, None)
            allmass.append(m)
            top = set(raw_data[m]["nodes"].keys())
            for chilren in raw_data[m]["edges"].values():
                top = top - set(chilren)
            top = list(top)
            self.alledges[m] = top

        self.allnodestype["00000001"] = "glycan"
        self.alledges["00000001"] = allmass

	self.allinedges = defaultdict(set)
	for pa,chs in self.alledges.items():
	    for ch in chs:
		self.allinedges[ch].add(pa)

        return raw_data

    def root(self):
        return "00000001"

    def nodes(self):
        for n in self.allnodestype.keys():
            yield n

    def parents(self, accession):
        for p in self.allinedges.get(accession, []):
            yield p

    def children(self, accession):
        for c in self.alledges.get(accession, []):
            yield c

    def level(self, accession):
	if accession in self.allnodestype:
	    return self.allnodestype[accession][1].lower()
        return None

    def get_molecularweight(self, accession):
	if accession in self.allnodestype:
	    return self.allnodestype[accession][0]
	return None

    def get_basecomposition(self, accession):
	if accession in self.allnodestype:
	    bcomp = self.allnodestype[accession][4]
	    if bcomp and bcomp in self.allnodestype:
	        return bcomp
	return None

    def get_composition(self, accession):
	if accession in self.allnodestype:
	    comp = self.allnodestype[accession][3]
	    if comp and comp in self.allnodestype:
	        return comp
	return None

    def get_topology(self, accession):
	if accession in self.allnodestype:
	    topo = self.allnodestype[accession][2]
	    if topo and topo in self.allnodestype:
	        return topo
	return None

    def get_iupac_composition(self, accession):
        return self.monosaccharide_count.get(accession, None)

    def get_iupac_composition_for_viewer(self, accession):
        return self.get_iupac_composition(accession)

    def get_iupac_composition_str_for_viewer(self, accession):
        s = ""
        d = self.get_iupac_composition_for_viewer(accession)
        for iupac in sorted(d.keys()):
            s += iupac+str(d[iupac])
        return s

    def get_missing_rank(self, accession):
        return self.allmissingrank[accession]

    def regexget(self, p, s):
        searchres = list(p.finditer(s))
        if len(searchres) == 1:
            return searchres[0].groupdict()["ans"]
        else:
            return None

    errorcategory = [
        ["cannot_parse_glycan_at_upper_subsumption",
         r"annotated (topology|composition|base composition) (?P<ans>(G\d{5}\w{2})) of G\d{5}\w{2} cannot be parsed"],
        ["cannot_parse_glycan_unknown", r"unknown problem parsing glycan (?P<ans>(G\d{5}\w{2}))"],
        ["unsupported_skeleton_code", r"unsupported skeleton code: (?P<ans>(.{1,10})) in glycan G\d{5}\w{2}"],
        ["mass_cannot_be_computed", r"mass could not be computed for (?P<ans>(G\d{5}\w{2}))"],
        ["mass_inconsistency",
         r"annotated mass \d*\.\d* for (?P<ans>(G\d{5}\w{2})) is different than computed mass \d*\.\d*"],
        # ["", r""],
    ]

    def geterrortype(self, s):
        for errortype, p in self.errorcategory:
            res0 = self.regexget(p, s)
            if res0:
                return errortype, res0
        return None, s

    def getwarningsbytype(self):
        if self.warningsbytype:
            return self.warningsbytype

        res = defaultdict(set)

        for i in reduce(lambda x, y: list(set(list(x) + list(y))), self.warnings.values()):
            errortype, ans = self.geterrortype(i)
            res[errortype].add(ans)

        self.warningsbytype = res
        return res

    def warningbytypetotal(self):
        warnings = self.getwarningsbytype()
        temp = map(lambda x: len(x), warnings.values())
        return reduce(lambda x, y: x + y, temp)

    def generateOWL(self, input_file_path, output_file_path, mass_lut_file_path, version=None, exact_sym=None, specific_sym=None):

        if specific_sym is None:
            specific_sym = {}
        if exact_sym is None:
            exact_sym = []

        all_syms = defaultdict(list)
        for e_sym in exact_sym:
            sym_type = "exact"
            for acc, sym0 in e_sym.items():
                all_syms[acc].append((sym0, sym_type))

        for sym_type, s_sym in specific_sym.items():
            for acc, sym0 in s_sym.items():
                all_syms[acc].append((sym0, sym_type))

        #for acc, syms in all_syms.items():
        #   print acc, syms

        self.loaddata(input_file_path)

        r = OWLWriter(mass_LUT_file_path=mass_lut_file_path, version=version)

        for mass in self.nodes():
            if not self.ismolecularweight(mass):
                continue
            mass = "%.2f" % float(mass)

            nodes = self.descendants(mass)

            r.addNode(mass, nodetype="molecularweight")

            for n in nodes:
                # molecular weight has a space between two words
                r.addNode(n, nodetype=self.level(n))
                node_obj = r.getNode(n)
                node_obj.set_basecomposition(self.get_basecomposition(n))
                node_obj.set_composition(self.get_composition(n))
                node_obj.set_topology(self.get_topology(n))
                node_obj.set_iupac_composition(self.get_iupac_composition_str_for_viewer(n))
                node_obj.set_missing_rank(self.get_missing_rank(n))


                if n in all_syms:
                    for sym0, sym0_type in all_syms[n]:
                        node_obj.addSynonym(sym0, sym0_type)

            nodes.add(mass)

            for n in nodes:
                if self.children(n):
                    for c in self.children(n):
                        r.connect(r.getNode(n), r.getNode(c), "subsumes")

        f = open(output_file_path, "w")
        s = r.write(f, r.make_graph())
        f.close()




from rdflib import URIRef, Namespace


class OWLWriter():
    _nodes = {}

    def __init__(self, mass_LUT_file_path=None, version=None):
        self.version = version
        if mass_LUT_file_path:
            self.mass_LUT_file_path = mass_LUT_file_path
        else:
            self.mass_LUT_file_path = "./mass_lut.txt"
        self.readmassidmap(mass_LUT_file_path)


    def addNode(self, nodeID, nodetype=None):
        self._nodes[nodeID] = NormalNode(nodeID, nodetype=nodetype)

    def addNodes(self, nodes):
        for nodeID in nodes:
            if nodeID not in self._nodes:
                self.addNode(nodeID)
            else:
                raise Exception

    def connect(self, n1, n2, edgeType):
        n1.connect(n2, edgeType)

    def getNode(self, nodeID):
        return self._nodes[nodeID]

    def allRelationship(self):
        res = set()
        for n in self._nodes.values():
            for r in n.getLink():
                res.add(r)
        res = list(res)
        return res

    def readmassidmap(self, mass_LUT_file_path):
        d = {}
        if mass_LUT_file_path:
            mass_lut_file_content = open(mass_LUT_file_path).read().strip().split("\n")
            for e, i in enumerate(mass_lut_file_content):
                if e == 0:
                    continue
                x = i.split("\t")
                d[float(x[1])] = x[0]
        else:
            d[1.01] = "10000001"
        self.massiddict = d

    newMass = False

    def overwritemasslookuptable(self):
        print "new mass ID was assigned"
        mass_lut_file_handle = open(self.mass_LUT_file_path, "w")
        mass_lut_file_handle.write("id\tmass\n")
        for mass in sorted(self.massiddict.keys()):
            id = self.massiddict[mass]
            mass_lut_file_handle.write("%s\t%.2f\n" % (id, mass))
        mass_lut_file_handle.close()

    gno = "http://purl.obolibrary.org/obo/"
    gtcs = "http://glytoucan.org/Structures/Glycans/"
    iao = "http://purl.obolibrary.org/obo/"
    gtco = "http://www.glytoucan.org/glyco/owl/glytoucan#"
    rocs = "http://www.glycoinfo.org/glyco/owl/relation#"
    oboInOwl = "http://www.geneontology.org/formats/oboInOwl#"

    glycan_class = 1
    definition_class = "IAO_0000115"
    subsumption_level_class = 11
    subsumption_level_annotation_property = 21
    glytoucan_id_annotation_property = 22
    glytoucan_link_annotation_property = 23

    structure_browser_link = 41
    composition_browser_link = 42

    cbbutton_annotation_property = 101

    has_Byonic_name_annotation_property = 202

    subsumption_level = {

        "molecularweight": {"id": 12,
                            "label": "molecular weight",
                            "definition": """
                                A subsumption category for glycans
                                described by their underivatized
                                molecular weight.
					  """,
                            "seeAlso": ["gtco:has_derivatization_type",
                                        "gtco:derivatization_type",
                                        "gtco:derivatization_type_permethylated"],
                            "comment": """
                                Underivatized molecular weight: The
                                molecular weight in the absence of any
                                chemical manipulation of a glycan for
                                analytical purposes. A common glycan
                                derivitization (chemical manipulation)
                                that affects glycan molecular weight
                                is permethylation.
                                       """
                            },

        "basecomposition": {"id": 13,
                            "label": "basecomposition",
                            "definition": """
				A subsumption category for glycans
				described by the number and type of
				monosaccharides with no monosaccharide
				stereochemistry or glycosidic bonds
				linking monosaccharides indicated.
					  """,
                            "seeAlso": ["rocs:Base_composition"]
                            },

        "composition": {"id": 14,
                        "label": "composition",
                        "definition": """
				A subsumption category for glycans
				described by the number and type
				of monosaccharides with partial or
				complete monosaccharide stereochemistry,
				but with no glycosidic bonds linking
				monosaccharides indicated.
					  """,
                        "seeAlso": ["rocs:Monosaccharide_composition"]
                        },

        "topology": {"id": 15,
                     "label": "topology",
                     "definition": """
				A subsumption category for glycans
				described by the arrangement of
				monosaccharides and the glycosidic bonds
				linking them, but with no linkage position
				or anomeric configuration indicated.
					  """,
                     "seeAlso": ["rocs:Glycosidic_topology"]
                     },

        "saccharide": {"id": 16,
                       "label": "saccharide",
                       "definition": """
                                A subsumption category for glycans
                                described by the arrangement of
                                monosaccharides and the glycosidic
                                bonds linking them, and with partial or
                                complete linkage position or anomeric
                                configuration indicated.
					  """,
                       "seeAlso": ["rocs:Linkage_defined_saccharide"]
                       },

    }

    has_subsumption_level = {

        "basecomposition": {
            "id": 33,
            "label": "has_basecomposition",
            "definition": """A metadata relation between a glycan and a glycan, with subsumption level basecomposition, that is equivalent to the glycan after the removal of all glycosidic bonds linking its monosaccharides and all monosaccharide stereochemistry information."""
        },

        "composition": {
            "id": 34,
            "label": "has_composition",
            "definition": """A metadata relation between a glycan and a glycan, with subsumption level composition, that is equivalent to the glycan after the removal of all glycosidic bonds linking its monosaccharides."""
        },

        "topology": {
            "id": 35,
            "label": "has_topology",
            "definition": """A metadata relation between a glycan and a glycan, with subsumption level topology, that is equivalent to the glycan after the removal of all linkage positions and anomeric configurations."""
        }
    }

    def gnoid(self, id):
        try:
            id = int(id)
            return "GNO_%08d" % (id,)
        except:
            pass
        assert re.search('^G[0-9]{5}[A-Z]{2}$', id)
        return "GNO_" + id

    def gnouri(self, id):
        return rdflib.URIRef(self.gno + self.gnoid(id))

    def make_graph(self):

        outputGraph = rdflib.Graph()

        rdf = rdflib.RDF
        rdfs = rdflib.RDFS
        owl = rdflib.OWL
        dcterms = rdflib.namespace.DCTERMS
        dc = rdflib.namespace.DC

        Literal = rdflib.Literal

        gtcs = Namespace(self.gtcs)
        iao = Namespace(self.iao)
        gno = Namespace(self.iao + "GNO_")
        gtco = Namespace(self.gtco)
        rocs = Namespace(self.rocs)
        oboInOwl = Namespace(self.oboInOwl)

        outputGraph.bind("owl", owl)
        outputGraph.bind("gno", gno)
        outputGraph.bind("obo", iao)
        outputGraph.bind("dc", dc)
        outputGraph.bind("dcterms", dcterms)
        outputGraph.bind("rocs", rocs)
        outputGraph.bind("oboInOwl", oboInOwl)

        root = rdflib.URIRef("http://purl.obolibrary.org/obo/gno.owl")

        # Add ontology
        outputGraph.add((root, rdf.type, owl.Ontology))
        outputGraph.add(
            (root, dc.title, Literal("Glycan Naming Ontology")))
        outputGraph.add(
            (root, dc.description, Literal("An ontology for glycans based on GlyTouCan, but organized by subsumption.")))

        # VersionIRI
        if self.version:
            outputGraph.add(
                (root, owl.versionIRI, URIRef("http://purl.obolibrary.org/obo/gno/%s/gno.owl" % str(datetime.date.today())))
            )

            outputGraph.add((root, owl.versionInfo, Literal("%s" % self.version)))

        # Copyright
        outputGraph.add(
            (root, dcterms.license, URIRef("http://creativecommons.org/licenses/by/4.0/")))
        outputGraph.add((root, rdfs.comment, Literal(
            "Glycan Naming Ontology is licensed under CC BY 4.0 (https://creativecommons.org/licenses/by/4.0/).")))
        outputGraph.add((root, rdfs.comment, Literal(" ".join("""
               Glycan Naming Ontology is licensed under CC BY 4.0. You
               are free to share (copy and redistribute the material
               in any medium or format) and adapt (remix, transform,
               and build upon the material) for any purpose, even
               commercially. for any purpose, even commercially. The
               licensor cannot revoke these freedoms as long as you follow
               the license terms. You must give appropriate credit (by
               using the original ontology IRI for the whole ontology and
               original term IRIs for individual terms), provide a link
               to the license, and indicate if any changes were made. You
               may do so in any reasonable manner, but not in any way
               that suggests the licensor endorses you or your use.
        """.split()))))

        # Add AnnotationProperty for definition
        definition = iao[self.definition_class]
        outputGraph.add((definition, rdf.type, owl.AnnotationProperty))
        outputGraph.add((definition, rdfs.isDefinedBy, iao["iao.owl"]))
        outputGraph.add((definition, rdfs.label, Literal("definition")))

        # Add AnnotationProperty for subsumption level
        has_subsumption_level_node = self.gnouri(self.subsumption_level_annotation_property)

        outputGraph.add((has_subsumption_level_node, rdf.type, owl.AnnotationProperty))
        outputGraph.add(
            (has_subsumption_level_node, rdfs.label, Literal("has_subsumption_category")))
        outputGraph.add((has_subsumption_level_node, definition,
                         Literal("A metadata relation between a glycan and its subsumption category.")))

        # Add AnnotationProperty for linking Glytoucan
        has_glytoucan_id_node = self.gnouri(self.glytoucan_id_annotation_property)

        outputGraph.add((has_glytoucan_id_node, rdf.type, owl.AnnotationProperty))
        outputGraph.add((has_glytoucan_id_node, rdfs.label, Literal("has_glytoucan_id")))
        outputGraph.add((has_glytoucan_id_node, definition,
                         Literal("The accession of the GlyTouCan entry describing the indicated glycan.")))

        has_glytoucan_link_node = self.gnouri(self.glytoucan_link_annotation_property)

        outputGraph.add((has_glytoucan_link_node, rdf.type, owl.AnnotationProperty))
        outputGraph.add((has_glytoucan_link_node, rdfs.label, Literal("has_glytoucan_link")))
        outputGraph.add((has_glytoucan_link_node, definition,
                         Literal("The URL of the GlyTouCan entry describing the indicated glycan.")))

        # Add sumbsumption level class and its instances
        rdfNodeSubsumption = self.gnouri(self.subsumption_level_class)

        outputGraph.add((rdfNodeSubsumption, rdf.type, owl.Class))
        outputGraph.add((rdfNodeSubsumption, rdfs.label, Literal("subsumption category")))
        outputGraph.add((rdfNodeSubsumption, definition,
                         Literal("Extent of glycan characterization provided by a glycan description.")))

        subsumptionLevel = {}

        for level in ("molecularweight", "basecomposition", "composition", "topology", "saccharide"):
            details = self.subsumption_level[level]
            rdfNode = self.gnouri(details["id"])
            subsumptionLevel[level] = rdfNode
            outputGraph.add((rdfNode, rdf.type, owl.NamedIndividual))
            outputGraph.add((rdfNode, rdf.type, rdfNodeSubsumption))
            outputGraph.add((rdfNode, rdfs.label, Literal(details["label"])))
            outputGraph.add((rdfNode, definition, Literal(" ".join(details["definition"].split()))))
            if "comment" in details:
                outputGraph.add((rdfNode, rdfs.comment, Literal(" ".join(details["comment"].split()))))
            for sa in details.get('seeAlso', []):
                ns, term = sa.split(":")
                outputGraph.add((rdfNode, rdfs.seeAlso, eval("%s.%s" % (ns, term))))

        # Add AnnotationProperty for 3 different has subsumption level
        has_xxx_nodes = {}
        for sl in ("basecomposition", "composition", "topology"):
            has_xxx_node = self.gnouri(self.has_subsumption_level[sl]["id"])

            outputGraph.add((has_xxx_node, rdf.type, owl.AnnotationProperty))
            outputGraph.add((has_xxx_node, rdfs.label, Literal(self.has_subsumption_level[sl]["label"])))
            outputGraph.add((has_xxx_node, definition, Literal(self.has_subsumption_level[sl]["definition"])))

            has_xxx_nodes[sl] = has_xxx_node

        has_topology_node = has_xxx_nodes["topology"]
        has_composition_node = has_xxx_nodes["composition"]
        has_basecomposition_node = has_xxx_nodes["basecomposition"]

        # Add AnnotationProperty for IUPAC composition (for the viewer)
        cbbutton_node = self.gnouri(self.cbbutton_annotation_property)

        outputGraph.add((cbbutton_node, rdf.type, owl.AnnotationProperty))
        outputGraph.add((cbbutton_node, rdfs.label, Literal("_widget_button_state")))

        # Add AnnotationProperty quick access to the GNOme browser
        has_structure_browser_node = self.gnouri(self.structure_browser_link)

        outputGraph.add((has_structure_browser_node, rdf.type, owl.AnnotationProperty))
        outputGraph.add((has_structure_browser_node, rdfs.label, Literal("has_structure_browser_link")))
        outputGraph.add((has_structure_browser_node, definition,
                         Literal("The URL of GNOme structure browser entry for the indicated glycan.")))

        has_composition_browser_node = self.gnouri(self.composition_browser_link)

        outputGraph.add((has_composition_browser_node, rdf.type, owl.AnnotationProperty))
        outputGraph.add((has_composition_browser_node, rdfs.label, Literal("has_composition_browser_link")))
        outputGraph.add((has_composition_browser_node, definition,
                         Literal("The URL of GNOme composition browser entry for the indicated glycan.")))



        # Add AnnotationProperty for hasExactSynonym
        hasExactSynonym_node = oboInOwl["hasExactSynonym"]

        outputGraph.add((hasExactSynonym_node, rdf.type, owl.AnnotationProperty))
        outputGraph.add((hasExactSynonym_node, rdfs.isDefinedBy, oboInOwl[""]))
        outputGraph.add((hasExactSynonym_node, rdfs.label, Literal("hasExactSynonym")))

        # Add AnnotationProperty for has_Byonic_name
        has_Byonic_name_node = self.gnouri(self.has_Byonic_name_annotation_property)

        outputGraph.add((has_Byonic_name_node, rdf.type, owl.AnnotationProperty))
        outputGraph.add((has_Byonic_name_node, definition, Literal("Glycan composition string used by the Byonic search engine.")))
        outputGraph.add((has_Byonic_name_node, rdfs.label, Literal("has_Byonic_name")))

        sym_types = {
            "exact": hasExactSynonym_node,
            "byonic": has_Byonic_name_node
        }

        # Glycan class under OWL things
        rdfNode = self.gnouri(self.glycan_class)

        outputGraph.add((rdfNode, rdf.type, owl.Class))
        outputGraph.add((rdfNode, rdfs.label, Literal("glycan")))
        outputGraph.add((rdfNode, definition,
                         Literal("A compound consisting of monosaccharides linked by glycosidic bonds.")))
        outputGraph.add((rdfNode, rdfs.seeAlso,
                         rdflib.URIRef("http://purl.obolibrary.org/obo/CHEBI_50699")))
        outputGraph.add((rdfNode, rdfs.seeAlso,
                         rdflib.URIRef("http://purl.obolibrary.org/obo/CHEBI_18154")))

        for n in self._nodes.values():
            if not n.ismolecularweight():
                rdfNode = self.gnouri(n.getID())
            else:
                try:
                    id = self.massiddict[float(n.getID())]
                except KeyError:
                    mass = n.getID()
                    id = str(int(max(self.massiddict.values())) + 1)
                    self.massiddict[float(mass)] = id
                    self.newMass = True
                rdfNode = self.gnouri(id)
            outputGraph.add((rdfNode, rdf.type, owl.Class))

            if n._nodeType:
                # outputGraph.add((rdfNode, rdfs.label, ns2[n._nodeType]))
                outputGraph.add((rdfNode, has_subsumption_level_node, subsumptionLevel[n._nodeType.lower()]))
            else:
                raise ValueError

            if not n.ismolecularweight():
                outputGraph.add((rdfNode, has_glytoucan_id_node, Literal(n.getID())))
                outputGraph.add((rdfNode, has_glytoucan_link_node, gtcs[n.getID()]))
                outputGraph.add((rdfNode, definition,
                                 Literal("A glycan described by the GlyTouCan entry with accession %s." % n.getID())))
                outputGraph.add((rdfNode, rdfs.label,
                                 Literal("%s" % n.getID())))
                # TODO add missing rank to GNOme latter here
                # print n.getID(), n.missing_rank()
            else:
                outputGraph.add((rdfNode, rdfs.subClassOf, self.gnouri(self.glycan_class)))
                outputGraph.add((rdfNode, rdfs.label,
                                 Literal("glycan of molecular weight %s Da" % n.getID())))
                outputGraph.add((rdfNode, definition,
                                 Literal(
                                     "A glycan characterized by underivitized molecular weight of %s Daltons." % n.getID())))

            for sym, sym_type in n.synonym():
                outputGraph.add((rdfNode, sym_types[sym_type], Literal(sym)))

            bcomp_acc = n.get_basecomposition()
            comp_acc = n.get_composition()
            topo_acc = n.get_topology()
            if bcomp_acc:
                outputGraph.add((rdfNode, has_basecomposition_node, self.gnouri(bcomp_acc)))
            if comp_acc:
                outputGraph.add((rdfNode, has_composition_node, self.gnouri(comp_acc)))
            if topo_acc:
                outputGraph.add((rdfNode, has_topology_node, self.gnouri(topo_acc)))

            cbbutton_str = n.get_iupac_composition()
            if cbbutton_str:
                outputGraph.add((rdfNode, cbbutton_node, Literal(cbbutton_str)))

            if not n.ismolecularweight():
                if n.isbasecomposition() or n.iscomposition():
                    outputGraph.add((rdfNode, has_composition_browser_node, URIRef(
                        "https://gnome.glyomics.org/GNOme.compositionselector.html?focus=%s" % n.getID())))
                outputGraph.add((rdfNode, has_structure_browser_node, URIRef(
                    "https://gnome.glyomics.org/GNOme.browser.html?focus=%s" % n.getID())))


        for l in self.allRelationship():
            if l.getEdgeType() == "subsumes":
                if l._sideA._nodeType == "molecularweight":
                    id = self.massiddict[float(l._sideA.getID())]
                    n1 = self.gnouri(id)
                else:
                    n1 = self.gnouri(l._sideA.getID())
                n2 = self.gnouri(l._sideB.getID())
                outputGraph.add((n2, rdfs.subClassOf, n1))

        if self.newMass:
            self.overwritemasslookuptable()

        return outputGraph

    def write(self, handle, graph):
        from rdflib.plugins.serializers.rdfxml import PrettyXMLSerializer
        writer = PrettyXMLSerializer(SubjectOrderedGraph(graph), max_depth=1)
        writer.serialize(handle)


# Defined solely to manipulate the accessor methods used by
# rdflib.plugins.serializers.rdfxml.PrettyXMLSerializer so that subjects
# and predictates are output in a deterministic order, and without nesting.

class SubjectOrderedGraph(object):

    def __init__(self, graph):
        self.graph = graph
        self.namespace_manager = self.graph.namespace_manager
        self.subjects_callno = 0
        self.subjects_order = {'gno.owl': 0, 'IAO_0000115': 1}

    def sortkey(self, *args):
        return tuple(map(str, args))

    def subject_sortkey(self, node):
        strnode = str(node)
        item = strnode.rsplit('/', 1)[1]
        return (self.subjects_order.get(item, 2), strnode)

    def subjects(self, *args, **kwargs):
        self.subjects_callno += 1
        if self.subjects_callno > 1:
            for n in sorted(self.graph.subjects(*args, **kwargs),
                            key=self.subject_sortkey):
                yield n

    def predicate_objects(self, *args, **kwargs):
        for p, o in sorted(self.graph.predicate_objects(*args, **kwargs),
                           key=lambda t: self.sortkey(*t)):
            yield p, o

    def objects(self, *args, **kwargs):
        for o in sorted(self.graph.objects(*args, **kwargs),
                        key=self.sortkey, reverse=True):
            yield o

    def __contains__(self, *args, **kwargs):
        return self.graph.__contains__(*args, **kwargs)

    def predicates(self, *args, **kwargs):
        return self.graph.predicates(*args, **kwargs)

    def triples_choices(self, *args, **kwargs):
        return self.graph.triples_choices(*args, **kwargs)



class NormalNode:
    _id = None
    _nodeType = None
    _nodeTypes = tuple(["saccharide", "topology", "composition", "basecomposition", "molecularweight"])
    _synonymTypes = tuple(["exact", "byonic"])


    def __init__(self, id, nodetype=None):
        self._links = []
        self._id = id
        self._nodeType = nodetype
        self._synonym = set()
        self._missing_rank = None

    def __str__(self):
        return self._id

    def getID(self):
        return self._id

    def getLink(self):
        return self._links

    def connect(self, node, edgeType=None):
        l = NormalLink(self, node, edgeType)

    def disconnect(self, node, edgeType=None):
        for l in self._links:
            if self in l.getBothNodes() and node in l.getBothNodes():  # and edgeType == l.getEdgeType:
                l.delete()

    def getparent(self):
        res = []
        for l in self._links:
            if l.getEdgeType() == "subsumes" and l._sideB == self:
                res.append(l.getoppositenode(self))
        return res

    def getchild(self):
        res = []
        for l in self._links:
            if l.getEdgeType() == "subsumes" and l._sideA == self:
                res.append(l.getoppositenode(self))
        return res

    def getequavalent(self):
        res = []
        for l in self._links:
            if l.getEdgeType() == "equal":
                res.append(l.getoppositenode(self))
        return res

    def getdirectlylinkednodes(self):
        res = []
        for l in self._links:
            res.append(l.getoppositenode(self))
        return res

    def set_basecomposition(self, bc):
        self._basecomposition = bc

    def get_basecomposition(self):
        try:
            return self._basecomposition
        except AttributeError:
            return None

    def set_composition(self, c):
        self._composition = c

    def get_composition(self):
        try:
            return self._composition
        except AttributeError:
            return None

    def set_topology(self, topology):
        self._topology = topology

    def get_topology(self):
        try:
            return self._topology
        except AttributeError:
            return None

    def set_iupac_composition(self, ic):
        self._iupac_composition_str = ic

    def get_iupac_composition(self):
        try:
            return self._iupac_composition_str
        except AttributeError:
            return None

    def addSynonym(self, sym, t):
        assert t in self._synonymTypes

        pair = (sym, t)
        self._synonym.add(pair)

        if t != "exact":
            pair = (sym, "exact")
            self._synonym.add(pair)

    def synonym(self):

        for sym in self._synonym:
            yield sym

    def ismolecularweight(self):
        return self._nodeType == "molecularweight"

    def isbasecomposition(self):
        return self._nodeType == "basecomposition"

    def iscomposition(self):
        return self._nodeType == "composition"

    def istopology(self):
        return self._nodeType == "topology"

    def issaccharide(self):
        return self._nodeType == "saccharide"

    def set_missing_rank(self, mr):
        self._missing_rank = mr

    def missing_rank(self):
        return self._missing_rank



class NormalLink:
    _id = None
    _edgeType = None
    _edgeTypes = tuple(["equal", "subsumes"])
    _exist = False
    # Real node object rather than node ID
    _sideA = None  # Parent
    _sideB = None  # Child

    def __init__(self, n1, n2, et):
        self._exist = True
        self._sideA = n1
        self._sideB = n2
        self._edgeType = et
        n1.getLink().append(self)
        n2.getLink().append(self)

    def __str__(self):
        return "%s -- %s --> %s" % (self._sideA.getID(), self._edgeType, self._sideB.getID())

    def getEdgeType(self):
        return self._edgeType

    def getBothNodes(self):
        return [self._sideA, self._sideB]

    def getoppositenode(self, nodeOBJ):
        if nodeOBJ == self._sideA:
            return self._sideB
        elif nodeOBJ == self._sideB:
            return self._sideA
        else:
            raise IndexError

    def delete(self):
        self._sideA.getLink().remove(self)
        self._sideB.getLink().remove(self)


# Theme generator for GNOme viewer
class GNOme_Theme_Base:

    #restriction_url = "https://raw.githubusercontent.com/glygen-glycan-data/GNOme/master/restrictions/GNOme_%s.accessions.txt"
    restriction_url = ""

    def __init__(self, restrictions_path):
        self.restriction_url = restrictions_path + "/GNOme_%s.accessions.txt"

    def getdata(self):
        # icon style
        # image source
        # external resources
        #   glycan URL
        #   Name
        #   Set
        raise NotImplemented()

    def get_accessions(self, restriction_set_name):
        # url = self.restriction_url%restriction_set_name
        # accs = urllib2.urlopen(url).read().strip().split()
        url = self.restriction_url % restriction_set_name
        accs = list(sorted(open(url).read().strip().split()))
        return accs

    def getoutputpath(self):
        raise NotImplemented()

    def write(self, output_path):
        d = self.getdata()
        json.dump(d, open(output_path, "w"))

class GNOme_Theme_GlyTouCan(GNOme_Theme_Base):

    def getdata(self):
        return {
            "icon_style": "snfg",
            "image_source_prefix": "https://image.glycosmos.org/snfg/png/",
            "image_source_suffix": "",
            "external_resources": [
                {
                    "name": "GlyTouCan",
                    "url_prefix": "https://glytoucan.org/Structures/Glycans/",
                    "url_suffix": "",
                    "glycan_set": None
                }
            ]

        }

class GNOme_Theme_GlyGen(GNOme_Theme_Base):

    def getdata(self):
        return {
            "icon_style": "snfg",
            "image_url_prefix": "https://api.glygen.org/glycan/image/",
            "image_url_suffix": "",
            "external_resources": [
                {
                    "name": "GlyGen",
                    "url_prefix": "https://www.glygen.org/glycan_detail.html?glytoucan_ac=",
                    "url_suffix": "",
                    "glycan_set": self.get_accessions("GlyGen")
                }
            ]

        }

class GNOme_Theme_GlyGenDev(GNOme_Theme_Base):

    def getdata(self):
        return {
            "icon_style": "snfg",
            "image_url_prefix": "https://api.glygen.org/glycan/image/",
            "image_url_suffix": "",
            "external_resources": [
                {
                    "name": "GlyGen Dev",
                    "url_prefix": "https://beta.glygen.org/glycan_detail.html?glytoucan_ac=",
                    "url_suffix": "",
                    "glycan_set": self.get_accessions("GlyGen")
                }
            ]

        }

class GNOme_Theme_Default(GNOme_Theme_Base):

    def getdata(self):
        return {
            "icon_style": "snfg",
            "image_url_prefix": "https://image.glycosmos.org/snfg/png/",
            "image_url_suffix": "",
            "external_resources": [
                {
                    "name": "GlycanData",
                    "url_prefix": "https://edwardslab.bmcb.georgetown.edu/glycandata/",
                    "url_suffix": "",
                    "glycan_set": self.get_accessions("GlycanData")
                },{
                    "name": "GlyGen",
                    "url_prefix": "https://www.glygen.org/glycan_detail.html?glytoucan_ac=",
                    "url_suffix": "",
                    "glycan_set": self.get_accessions("GlyGen")
                },{
                    "name": "GlyTouCan",
                    "url_prefix": "https://glytoucan.org/Structures/Glycans/",
                    "url_suffix": "",
                    "glycan_set": None
                }
            ]

        }

import pygly.Monosaccharide
class IncompleteScore:

    def __init__(self):
        pass

    def monoscore(self, m):
        res = 0.0
        anomer = m.anomer()
        config = m.config()
        rs = m.ring_start()
        re = m.ring_end()
        stem = m.stem()
        superclass = m.superclass()

        # sublink = m.substituent_links()

        if anomer == None:
            res += 2

        if config != None:
            res += float(config.count(None)) / len(config) * 2
        else:
            res += 2

        if anomer != pygly.Monosaccharide.Anomer.uncyclized:
            if rs == None:
                res += 1
            if re == None:
                res += 1

        if stem == None or None in stem:
            res += 2

        if superclass == None:
            res += 2

        res = res # max 10
        return res

    def linksimplescore(self, l):
        res = 0.0
        pp = l.parent_pos()
        cp = l.child_pos()
        pt = l.parent_type()
        ct = l.child_type()
        und = l.undetermined()

        if pp != None:
            if len(pp) > 1:
                res += 1
        else:
            res += 2

        if cp != None:
            if len(cp) > 1:
                res += 1
        else:
            res += 2

        if len(pt) > 1:
            res += 2

        if len(ct) > 1:
            res += 2
        # print pp, cp, pt, ct, und

        res = res / 8 * 10 # scale from 8 to 10
        return res


    def substscore(self, s):
        if s._sub != None:
            return 2
        return 0

    def score(self, g):

        monos = list(g.all_nodes())
        total_mono = len(monos)

        unconnected_mono_root = list(g.unconnected_roots())
        unconnected_mono_root = filter(lambda m: m.is_monosaccharide(), unconnected_mono_root)
        det_parent_links = [m.parent_links() for m in filter(lambda m: m not in unconnected_mono_root and m != g.root(), monos)]
        und_parent_links = [m.parent_links() for m in unconnected_mono_root]

        allsubst = filter(lambda n: not n.is_monosaccharide(), g.all_nodes(subst=True, undet_subst=False))

        mono_score_result = 0.0
        mono_score_total  = 10.0 * total_mono

        link_score_result = 0.0
        link_score_total  = 10.0 * (total_mono - 1)

        for m in monos:
            mono_score_result += self.monoscore(m)

        if len(und_parent_links) == total_mono:
            # Basecomp...
            link_score_result = 10.0 * (total_mono - 1)
        else:
            for ls in und_parent_links + det_parent_links:

                if len(ls) == 1:
                    link_score_result += self.linksimplescore(ls[0]) * 0.25

                elif len(ls) > 1:
                    parent_link_scores = []
                    for pl in ls:
                        s = self.linksimplescore(pl) * 0.25
                        parent_link_scores.append(s)

                    ave = sum(parent_link_scores) / len(parent_link_scores)
                    und_ratio = (float(len(parent_link_scores)-1) / (total_mono))
                    s = ave + 7.5 * und_ratio
                    link_score_result += s

                else:
                    # Unknown of from-and-to and detail of the link
                    link_score_result += 10.0




        for subst in allsubst:
            if str(subst) == "anhydro":
                continue

            pl = subst.parent_links()
            link_score_total += 5

            if len(pl) == 1:
                link_score_result += self.linksimplescore(pl[0]) * 0.125

            elif len(pl) > 1:

                parent_link_scores = []
                for pl0 in pl:
                    s = self.linksimplescore(pl0) * 0.125
                    parent_link_scores.append(s)

                ave = sum(parent_link_scores) / len(parent_link_scores)
                und_ratio = (float(len(parent_link_scores)-1) / (total_mono))
                s = ave + 3.75 * und_ratio
                link_score_result += s

            else:
                link_score_result += 5.0


        if link_score_total ==0:
            if link_score_result != 0:
                raise RuntimeError
            link_score_total = 1

        res = mono_score_result/mono_score_total * 0.35 + link_score_result/link_score_total * 0.65
        #print mono_score_result / mono_score_total, link_score_result / link_score_total
        res = int(res * 10000)
        return res


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

        g = SubsumptionGraph()
        g.compute(*sys.argv[1:], verbose=verbose)

    elif cmd == "writeowl":
        # python GNOme.py writeowl ../smw/glycandata/data/gnome_subsumption_raw.txt ./GNOme.owl mass_lookup_2decimal v1.1.5

        # "../smw/glycandata/data/gnome_subsumption_raw.txt"

        kv_para = {
            "version": None
        }
        if len(sys.argv) < 4:
            print "Please provide dumpfile, output file path(with file name), mass LUT path and version (optional)"
            sys.exit(1)
        if len(sys.argv) > 4:
            versionTag = sys.argv[4]
            for k, v in zip(sys.argv[4::2], sys.argv[5::2]):#sys.argv[4::2]:
                assert k.startswith("-")
                k = k[1:]
                kv_para[k] = v

        def symFile2dict(fpath):
            res = {}
            tmp = open(fpath).read().strip().split("\n")
            for pair in tmp:
                k,v = pair.split()
                res[k] = v
            return res

        exactSym = []
        specificSym = {}

        for k, v in kv_para.items():
            if "_sym" not in k:
                continue

            k = k.split("_sym")[0]
            if k == "exact":
                exactSym.append(symFile2dict(v))
            else:
                specificSym[k] = symFile2dict(v)


        #for syms in exactSym:
        #    print len(syms)
        #for k, syms in specificSym.items():
        #    print k, len(syms)

        ifn = sys.argv[1]  # path to input file
        ofn = sys.argv[2]
        mass_lut = sys.argv[3]

        subsumption_instance = SubsumptionGraph()
        subsumption_instance.generateOWL(ifn, ofn, mass_lut, version=kv_para["version"], exact_sym=exactSym, specific_sym=specificSym)

        if "allExactSymOutput" in kv_para:
            allexactsym = {}
            for symset in exactSym + specificSym.values():
                for acc, sym0 in symset.items():
                    if acc not in allexactsym:
                        allexactsym[acc] = []
                    allexactsym[acc].append(sym0)

            allExactSymOutputF = open(kv_para["allExactSymOutput"], "w")
            for acc in sorted(allexactsym.keys()):
                for sym0 in allexactsym[acc]:
                    allExactSymOutputF.write("%s\t%s\n" % (acc, sym0))
            allExactSymOutputF.close()


    elif cmd == "writeresowl":

        if len(sys.argv) < 4:
            print "Please provide GNOme.owl, restriction set name, output file path"
            sys.exit(1)

        ifn = sys.argv[1]
        ofn = sys.argv[3]
        restriction_accs_file = sys.argv[2]

        accs = open(restriction_accs_file).read().strip().split()

        GNOme_res = GNOme(resource=ifn)
        GNOme_res.restrict(accs)

        f = open(ofn, "w")
        GNOme_res.write(f)
        f.close()

    elif cmd == "viewerdata":
        # python GNOme.py viewerdata ./GNOme.owl ./gnome_subsumption_raw.txt ./GNOme.browser.js

        if len(sys.argv) < 3:
            print "Please provide GNOme.owl and output file path"
            sys.exit(1)

        ifn = sys.argv[1]
        ofn_data = sys.argv[2]
        gnome = GNOme(resource=ifn)
        gnome.toViewerData(ofn_data)


    elif cmd == "UpdateAcc":

        if len(sys.argv) < 4:
            print "Please provide restriction set name and file path"
            sys.exit(1)

        restriction_set_name = sys.argv[1]
        fp = sys.argv[2]

        restriction_set = []
        if restriction_set_name == "BCSDB":
            sys.exit(0)

        #elif restriction_set_name == "GlyGen":
        #    fp = os.path.dirname(os.path.abspath(__file__)) + "/../smw/glycandata/data/glygen_accessions.txt"
        #    restriction_set = open(fp).read().strip().split()

        elif restriction_set_name in ["GlycanData", "GlyGen"]:
            glycandata_tsv_fp = "../smw/glycandata/export/allglycan.tsv"

            restriction_set = open(glycandata_tsv_fp).read().strip().split()
            restriction_set.pop(0)
        else:
            print "Restriction set: %s is not supported"
            sys.exit(1)

        open(fp, "w").write("\n".join(restriction_set))

        json_fp = open(sys.argv[3], "w")
        json.dump(restriction_set, json_fp, sort_keys=True, indent=2)

    elif cmd == "UpdateTheme":

        if len(sys.argv) < 3:
            print "Please provide restriction set name and file path"
            sys.exit(1)

        restriction_url = sys.argv[1]
        theme_path = sys.argv[2]

        td = GNOme_Theme_Default(restriction_url)
        tgtc = GNOme_Theme_GlyTouCan(restriction_url)
        tgg = GNOme_Theme_GlyGen(restriction_url)
        tggd = GNOme_Theme_GlyGenDev(restriction_url)

        td.write(theme_path + "default.json")
        tgtc.write(theme_path + "GlyTouCan.json")
        tgg.write(theme_path + "GlyGen.json")
        tggd.write(theme_path + "GlyGenDev.json")


    else:

	gnome = GNOme()

	if hasattr(gnome,cmd):
	    result = getattr(gnome,cmd)(*sys.argv[1:])
	    if isinstance(result,basestring):
		print result
	    else:
	        for r in result:
		    print r
	else:
            print >> sys.stderr, "Bad command: %s" % (cmd,)
            sys.exit(1)
