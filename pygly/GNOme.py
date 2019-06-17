import sys, os, os.path, urllib, rdflib, copy, re
from collections import defaultdict


class GNOme(object):
    version = "1.0.1"
    referenceowl = "https://raw.githubusercontent.com/glygen-glycan-data/GNOme/V%s/GNOme.owl" % (version,)
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
        self.ns['gno'] = rdflib.Namespace('http://ontology.glygen.org/gnome/GNO_')
        self.ns['rdf'] = rdflib.Namespace('http://www.w3.org/1999/02/22-rdf-syntax-ns#')
        self.ns['rdfs'] = rdflib.Namespace('http://www.w3.org/2000/01/rdf-schema#')
        self.ns[None] = rdflib.Namespace("")

    def triples(self, subj=None, pred=None, obj=None):
        for (s, p, o) in self.gnome.triples(((self.uri(subj) if subj else None), (self.uri(pred) if pred else None),
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
            else:
                olab = self.accession(o)
            attr[plab] = olab
        attr[u'level'] = self.level(accession)
        return attr

    def edges(self):
        for n in self.nodes():
            for p in self.parents(n):
                yield p, n

    def root(self):
        return "00000001"

    def parents(self, accession):
        uri = "gno:%s" % (accession,)
        for s, p, o in self.triples(uri, 'rdfs:subClassOf'):
            yield self.accession(o)

    def ancestors(self, accession):
        anc = set()
        for p in self.parents(accession):
            anc.add(p)
            anc.update(self.ancestors(p))
        return anc

    def children(self, accession):
        uri = "gno:%s" % (accession,)
        for s, p, o in self.triples(None, 'rdfs:subClassOf', uri):
            yield self.accession(s)

    def descendants(self, accession):
        desc = set()
        for c in self.children(accession):
            desc.add(c)
            desc.update(self.descendants(c))
        return desc

    def isleaf(self, accession):
        for ch in self.children(accession):
            return False
        return True

    def isroot(self, accession):
        return accession == self.root()

    def level(self, accession):
        uri = "gno:%s" % (accession,)
        for s, p, o in self.triples(uri, "gno:00000021", None):
            return " ".join(unicode(self.label(o)).split()[2:])

    def islevel(self, accession, level):
        return self.level(accession) == level

    def ismolecularweight(self, accession):
        return self.islevel(accession, 'molecular weight')

    def isbasecomposition(self, accession):
        return self.islevel(accession, 'basecomposition')

    def iscomposition(self, accession):
        return self.islevel(accession, 'composition')

    def istopology(self, accession):
        return self.islevel(accession, 'topology')

    def issaccharide(self, accession):
        return self.islevel(accession, 'saccharide')

    def get_basecomposition(self, accession):
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

    def get_composition(self, accession):
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

    def get_topology(self, accession):
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

    def has_basecomposition(self, accession):
        assert self.isbasecomposition(accession)
        return self.descendants(accession)

    def has_composition(self, accession):
        assert self.iscomposition(accession)
        return self.descendants(accession)

    def has_topology(self, accession):
        assert self.istopology(accession)
        return self.descendants(accession)

    def restrict(self, restriction):
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
                        toremove.add((n1, n3))
        for n1, n3 in toremove:
            parents[n1].remove(n3)

        for n in self.nodes():
            if n not in keep:
                self.gnome.remove((self.uri("gno:" + n), None, None))
            else:
                self.gnome.remove((self.uri("gno:" + n), self.uri("rdfs:subClassOf"), None))
                for n1 in parents[n]:
                    self.gnome.add((self.uri("gno:" + n), self.uri("rdfs:subClassOf"), self.uri("gno:" + n1)))

    def write(self, handle):
        self.gnome.serialize(handle)

    def dump(self):
        for acc in sorted(g.nodes()):
            print acc
            for k, v in sorted(g.attributes(acc).items()):
                print "  %s: %s" % (k, v)


from alignment import GlycanSubsumption, GlycanEqual
from GlyTouCan import GlyTouCan
from Monosaccharide import Anomer
import time


class SubsumptionGraph:
    def __init__(self, *args, **kwargs):
        self.gtc = GlyTouCan(usecache=True)
        self.subsumption = GlycanSubsumption()
        self.geq = GlycanEqual()
        self.verbose = kwargs.get('verbose', 0)

        masscluster = defaultdict(dict)
        if len(args) > 0:
            for mass in map(float, args):
                rmass = round(mass, 2)
                for glyacc in self.gtc.hasmass(rmass):
                    masscluster[rmass][glyacc] = dict(accession=glyacc)
        else:
            for glyacc, mass in self.gtc.allmass():
                rmass = round(mass, 2)
                masscluster[rmass][glyacc] = dict(accession=glyacc)
        for rmass, cluster in sorted(masscluster.items()):
            self.compute_component(rmass, cluster)

    def warning(self, msg, level):
        if self.verbose >= level:
            print "# WARNING:%d - %s" % (level, msg)

    def compute_component(self, rmass, cluster):
        start = time.time()

        print "# START %s - %d accessions in molecular weight cluster for %s" % (time.ctime(), len(cluster), rmass)
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
                        self.warning("unsupported skeleton code: " + skel + " in glycan " + acc, 2)
                else:
                    self.warning("unknown problem parsing glycan " + acc, 2)
                continue
            cluster[acc]['glycan'] = gly
            cluster[acc]['mass'] = self.gtc.getmass(acc)

        clusteracc = set(map(lambda t: t[0], filter(lambda t: t[1].has_key('glycan'), cluster.items())))

        outedges = defaultdict(set)
        inedges = defaultdict(set)

        for acc1 in sorted(clusteracc):
            gly1 = cluster[acc1]['glycan']
            for acc2 in sorted(clusteracc):
                gly2 = cluster[acc2]['glycan']
                if acc1 != acc2:
                    if self.subsumption.leq(gly1, gly2):
                        if not self.geq.eq(gly1, gly2) or acc2 < acc1:
                            outedges[acc2].add(acc1)
                            inedges[acc1].add(acc2)

        for acc in clusteracc:
            topo = self.gtc.gettopo(acc)
            if topo:
                if topo not in cluster:
                    self.warning("annotated topology %s of %s is not in %s rounded mass cluster" % (topo, acc, rmass),
                                 1)
                elif not cluster[topo].get('glycan'):
                    self.warning("annotated topology %s of %s cannot be parsed" % (topo, acc), 1)
                elif acc not in outedges[topo] and acc != topo:
                    self.warning("annotated topology %s does not subsume %s" % (topo, acc), 1)
            comp = self.gtc.getcomp(acc)
            if comp:
                if comp not in cluster:
                    self.warning(
                        "annotated composition %s of %s is not in %s rounded mass cluster" % (comp, acc, rmass), 1)
                elif not cluster[comp].get('glycan'):
                    self.warning("annotated composition %s of %s cannot be parsed" % (comp, acc), 1)
                elif acc not in outedges[comp] and acc != comp:
                    self.warning("annotated composition %s does not subsume %s" % (comp, acc), 1)
            bcomp = self.gtc.getbasecomp(acc)
            if bcomp:
                if bcomp not in cluster:
                    self.warning(
                        "annotated base composition %s of %s is not in %s rounded mass cluster" % (bcomp, acc, rmass),
                        1)
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
            elif abs(cluster[acc]['mass'] - umw) > 0.0001:
                self.warning(
                    "annotated mass %s for %s is different than computed mass %s" % (cluster[acc]['mass'], acc, umw), 1)

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
                    self.warning(
                        "annotation inferred level %s for %s != computed level Saccharide (anomer)" % (level, acc), 1)
                continue
            if self.any_parent_pos(gly):
                if not level:
                    cluster[acc]['level'] = 'Saccharide*'
                elif level != "Saccharide":
                    self.warning(
                        "annotation inferred level %s for %s != computed level Saccharide (parent_pos)" % (level, acc),
                        1)
                continue
            if self.any_links(gly):
                if not level:
                    cluster[acc]['level'] = 'Topology*'
                elif level != "Topology":
                    self.warning("annotation inferred level %s for %s != computed level Topology" % (level, acc), 1)
                continue
            if self.monosaccharide_count(gly) == 1 and self.any_ring(gly):
                if not level:
                    cluster[acc]['level'] = 'Topology*'
                elif level != "Topology":
                    self.warning("annotation inferred level %s for %s != computed level Topology" % (level, acc), 1)
                continue
            if self.any_stem(gly):
                if not level:
                    cluster[acc]['level'] = 'Composition*'
                elif level != "Composition":
                    self.warning("annotation inferred level %s for %s != computed level Composition" % (level, acc), 1)
                continue
            if not level:
                cluster[acc]['level'] = 'BaseComposition*'
            elif level != "BaseComposition":
                self.warning("annotation inferred level %s for %s != computed level BaseComposition" % (level, acc), 1)

        for acc in clusteracc:
            g = cluster[acc]
            if not g.get('topo'):
                if g.get('level') in ("Topology", "Topology*"):
                    g['topo'] = acc + "*"
                else:
                    for acc1 in inedges[acc]:
                        if cluster[acc1].get('level') in ("Topology", "Topology*"):
                            if not g.get('topo') or g.get('topo').rstrip("*") in inedges[acc1]:
                                g['topo'] = acc1 + "*"
            if not g.get('comp'):
                if g.get('level') in ("Composition", "Composition*"):
                    g['comp'] = acc + "*"
                else:
                    for acc1 in inedges[acc]:
                        if cluster[acc1].get('level') in ("Composition", "Composition*"):
                            if not g.get('comp') or g.get('comp').rstrip("*") in inedges[acc1]:
                                g['comp'] = acc1 + "*"
            if not g.get('bcomp'):
                if g.get('level') in ("BaseComposition", "BaseComposition*"):
                    g['bcomp'] = acc + "*"
                else:
                    for acc1 in inedges[acc]:
                        if cluster[acc1].get('level') in ("BaseComposition", "BaseComposition*"):
                            if not g.get('bcomp') or g.get('bcomp').rstrip("*") in inedges[acc1]:
                                g['bcomp'] = acc1 + "*"
        if len(clusteracc) < 1:
            print "# DONE - Elapsed time %.2f sec." % (time.time() - start,)
            sys.stdout.flush()
            return

        print "# NODES - %d/%d glycans in molecular weight cluster for %s" % (len(clusteracc), total, rmass)
        for acc in sorted(clusteracc, key=lambda acc: cluster[acc].get('level')):
            g = cluster[acc]
            print acc, g.get('mass'), g.get('level'), g.get('topo'), g.get('comp'), g.get('bcomp'),
            gly = g.get('glycan')
            print ("COMP " if (not gly.has_root()) else "") + \
                  ("UNDET " if (gly.undetermined() and gly.has_root()) else "") + \
                  ("FULL" if gly.fully_determined() else "")
        print "# ENDNODES - %d/%d glycans in molecular weight cluster for %s" % (len(clusteracc), total, rmass)
        sys.stdout.flush()

        prunedoutedges = self.prune_edges(outedges)
        prunedinedges = self.prune_edges(inedges)

        print "# TREE"
        for r in clusteracc:
            if len(prunedinedges[r]) == 0:
                self.print_tree(cluster, prunedoutedges, r)
        print "# ENDTREE"

        print "# EDGES"
        for n in clusteracc:
            print "%s:" % (n,),
            print " ".join(prunedoutedges[n])
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
            if l.parent_pos() != None and l.parent_pos() != set([1]):
                return True
        return False

    def any_ring(self, gly):
        for m in gly.all_nodes():
            if m.ring_start() or m.ring_end():
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

    def monosaccharide_count(self, gly):
        return sum(1 for _ in gly.all_nodes())

    def prune_edges(self, inedges):
        edges = copy.deepcopy(inedges)
        toremove = set()
        for n1 in list(edges):
            for n2 in edges[n1]:
                for n3 in edges[n2]:
                    assert n3 in edges[n1], ", ".join([n1, n2, n3, ":"] + edges[n1] + [":"] + edges[n2])
                    if n3 in edges[n1]:
                        toremove.add((n1, n3))
        for n1, n3 in toremove:
            edges[n1].remove(n3)
        return edges

    def print_tree(self, cluster, edges, root, indent=0):
        print "%s%s" % (" " * indent, root), cluster[root]['level']
        for ch in edges[root]:
            self.print_tree(cluster, edges, ch, indent + 2)


class NormalNode:
    _id = None
    _nodeType = None
    _nodeTypes = tuple(["saccharide", "topology", "composition", "basecomposition", "molecularweight"])

    def __init__(self, id, nodetype=None):
        self._links = []
        self._id = id
        self._nodeType = nodetype

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


class relationship():
    _nodes = {}

    def __init__(self):
        self.readmassidmap()

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

    def readmassidmap(self):
        d = {}
        # todo get from internet
        f = open("mass_lookup_2decimal").read().strip().split("\n")
        for e, i in enumerate(f):
            if e == 0:
                continue
            x = i.split("\t")
            d[float(x[1])] = x[0]
        self.massiddict = d

    newMass = False

    def overwritemasslookuptable(self):
        f = open("mass_lookup_2decimal_new", "w")
        for mass in sorted(self.massiddict.keys()):
            id = self.massiddict[mass]
            f.write("%s\t%.2f\n" % (id, mass))
        f.close()

    def OWLgenerate(self):
        outputGraph = rdflib.Graph()

        ns = rdflib.Namespace("http://ontology.glygen.org/gnome/")
        nsgtc = rdflib.Namespace("http://glytoucan.org/Structures/Glycans/")
        nsiao = rdflib.Namespace("http://purl.obolibrary.org/obo/")
        nsgtco = rdflib.Namespace("http://www.glytoucan.org/glyco/owl/glytoucan#")
        # nsrocks = rdflib.Namespace("")

        outputGraph.bind("owl", rdflib.OWL)
        outputGraph.bind("gno", ns)
        outputGraph.bind("obo", nsiao)
        outputGraph.bind("dc", rdflib.namespace.DC)
        # outputGraph.bind("rocks", nsrocks)

        root = rdflib.URIRef("http://purl.obolibrary.org/obo/gno.owl")

        # Add ontology
        outputGraph.add((root, rdflib.namespace.RDF.type, rdflib.OWL.Ontology))

        # Copyright
        outputGraph.add(
            (root, rdflib.namespace.DC.license, rdflib.URIRef("http://creativecommons.org/licenses/by/4.0/")))
        outputGraph.add((root, rdflib.namespace.RDFS.comment, rdflib.Literal(
            "Glycan Naming Ontology is licensed under CC BY 4.0 (https://creativecommons.org/licenses/by/4.0/).")))
        outputGraph.add((root, rdflib.namespace.RDFS.comment, rdflib.Literal(
            "Glycan Naming Ontology is licensed under CC BY 4.0. You are free to share (copy and redistribute the material in any medium or format) and adapt (remix, transform, and build upon the material) for any purpose, even commercially. for any purpose, even commercially. The licensor cannot revoke these freedoms as long as you follow the license terms. You must give appropriate credit (by using the original ontology IRI for the whole ontology and original term IRIs for individual terms), provide a link to the license, and indicate if any changes were made. You may do so in any reasonable manner, but not in any way that suggests the licensor endorses you or your use.")))

        # Add AnnotationProperty for definition
        definition_node = nsiao["IAO_0000115"]
        outputGraph.add((definition_node, rdflib.namespace.RDF.type, rdflib.OWL.AnnotationProperty))
        outputGraph.add((definition_node, rdflib.namespace.RDFS.isDefinedBy, nsiao["iao.owl"]))
        outputGraph.add((definition_node, rdflib.namespace.RDFS.label, rdflib.Literal("definition")))

        # Add AnnotationProperty for subsumption level
        has_subsumption_level_node = ns["GNO_00000021"]
        outputGraph.add((has_subsumption_level_node, rdflib.namespace.RDF.type, rdflib.OWL.AnnotationProperty))
        outputGraph.add(
            (has_subsumption_level_node, rdflib.namespace.RDFS.label, rdflib.Literal("has_subsumption_category")))
        outputGraph.add((has_subsumption_level_node, definition_node,
                         rdflib.Literal("A metadata relation between a glycan and its subsumption category.")))

        # Add AnnotationProperty for linking Glytoucan
        has_glytoucan_id_node = ns["GNO_00000022"]
        outputGraph.add((has_glytoucan_id_node, rdflib.namespace.RDF.type, rdflib.OWL.AnnotationProperty))
        outputGraph.add((has_glytoucan_id_node, rdflib.namespace.RDFS.label, rdflib.Literal("has_glytoucan_id")))
        outputGraph.add((has_glytoucan_id_node, definition_node,
                         rdflib.Literal("The accession of the GlyTouCan entry describing the indicated glycan.")))

        has_glytoucan_link_node = ns["GNO_00000023"]
        outputGraph.add((has_glytoucan_link_node, rdflib.namespace.RDF.type, rdflib.OWL.AnnotationProperty))
        outputGraph.add((has_glytoucan_link_node, rdflib.namespace.RDFS.label, rdflib.Literal("has_glytoucan_link")))
        outputGraph.add((has_glytoucan_link_node, definition_node,
                         rdflib.Literal("The URL of the GlyTouCan entry describing the indicated glycan.")))

        # Add sumbsumption level class and its instances
        rdfNodeSubsumption = ns["GNO_00000011"]
        outputGraph.add((rdfNodeSubsumption, rdflib.namespace.RDF.type, rdflib.OWL.Class))
        outputGraph.add((rdfNodeSubsumption, rdflib.namespace.RDFS.label, rdflib.Literal("subsumption category")))
        outputGraph.add((rdfNodeSubsumption, definition_node,
                         rdflib.Literal("Extent of glycan characterization provided by a glycan description.")))

        subsumptionLevel = {}
        rdfNode = ns["GNO_000000%s" % (12)]
        subsumptionLevel["molecularweight"] = rdfNode
        outputGraph.add((rdfNode, rdflib.namespace.RDF.type, rdflib.OWL.NamedIndividual))
        outputGraph.add((rdfNode, rdflib.RDF.type, rdfNodeSubsumption))
        outputGraph.add(
            (rdfNode, rdflib.namespace.RDFS.label, rdflib.Literal("subsumption category %s" % "molecular weight")))
        outputGraph.add((rdfNode, definition_node, rdflib.Literal(
            "A subsumption category for glycans described by their underivatized molecular weight.")))
        outputGraph.add((rdfNode, rdflib.namespace.RDFS.comment, rdflib.Literal(
            "Underivatized molecular weight: The molecular weight in the absence of any chemical manipulation of a glycan for analytical purposes. A common glycan derivitization (chemical manipulation) that affects glycan molecular weight is permethylation.")))
        outputGraph.add((rdfNode, rdflib.namespace.RDFS.seeAlso, nsgtco["has_derivatization_type"]))
        outputGraph.add((rdfNode, rdflib.namespace.RDFS.seeAlso, nsgtco["derivatization_type"]))
        # outputGraph.add((rdfNode, rdflib.namespace.RDFS.seeAlso, nsgtco["derivatization_type_none"]))
        outputGraph.add((rdfNode, rdflib.namespace.RDFS.seeAlso, nsgtco["derivatization_type_permethylated"]))
        # outputGraph.add((rdfNode, rdflib.namespace.RDFS.seeAlso, nsgtco["derivatization_type_peracetylated"]))

        rdfNode = ns["GNO_000000%s" % (13)]
        subsumptionLevel["basecomposition"] = rdfNode
        outputGraph.add((rdfNode, rdflib.namespace.RDF.type, rdflib.OWL.NamedIndividual))
        outputGraph.add((rdfNode, rdflib.RDF.type, rdfNodeSubsumption))
        outputGraph.add(
            (rdfNode, rdflib.namespace.RDFS.label, rdflib.Literal("subsumption category %s" % "basecomposition")))
        outputGraph.add((rdfNode, definition_node, rdflib.Literal(
            "A subsumption category for glycans described by the number and type of monosaccharides with no monosaccharide stereochemistry or glycosidic bonds linking monosaccharides indicated.")))
        outputGraph.add((rdfNode, rdflib.namespace.RDFS.seeAlso,
                         rdflib.URIRef("http://www.glycoinfo.org/glyco/owl/relation#Base_composition")))

        rdfNode = ns["GNO_000000%s" % (14)]
        subsumptionLevel["composition"] = rdfNode
        outputGraph.add((rdfNode, rdflib.namespace.RDF.type, rdflib.OWL.NamedIndividual))
        outputGraph.add((rdfNode, rdflib.RDF.type, rdfNodeSubsumption))
        outputGraph.add(
            (rdfNode, rdflib.namespace.RDFS.label, rdflib.Literal("subsumption category %s" % "composition")))
        outputGraph.add((rdfNode, definition_node, rdflib.Literal(
            "A subsumption category for glycans described by the number and type of monosaccharides with partial or complete monosaccharide stereochemistry, but with no glycosidic bonds linking monosaccharides indicated.")))
        outputGraph.add((rdfNode, rdflib.namespace.RDFS.seeAlso,
                         rdflib.URIRef("http://www.glycoinfo.org/glyco/owl/relation#Monosaccharide_composition")))

        rdfNode = ns["GNO_000000%s" % (15)]
        subsumptionLevel["topology"] = rdfNode
        outputGraph.add((rdfNode, rdflib.namespace.RDF.type, rdflib.OWL.NamedIndividual))
        outputGraph.add((rdfNode, rdflib.RDF.type, rdfNodeSubsumption))
        outputGraph.add((rdfNode, rdflib.namespace.RDFS.label, rdflib.Literal("subsumption category %s" % "topology")))
        outputGraph.add((rdfNode, definition_node, rdflib.Literal(
            "A subsumption category for glycans described by the arrangement of monosaccharides and the glycosidic bonds linking them, but with no linkage position or anomeric configuration indicated.")))
        outputGraph.add((rdfNode, rdflib.namespace.RDFS.seeAlso,
                         rdflib.URIRef("http://www.glycoinfo.org/glyco/owl/relation#Glycosidic_topology")))

        rdfNode = ns["GNO_000000%s" % (16)]
        subsumptionLevel["saccharide"] = rdfNode
        outputGraph.add((rdfNode, rdflib.namespace.RDF.type, rdflib.OWL.NamedIndividual))
        outputGraph.add((rdfNode, rdflib.RDF.type, rdfNodeSubsumption))
        outputGraph.add((rdfNode, rdflib.namespace.RDFS.label, rdflib.Literal("subsumption level %s" % "saccharide")))
        outputGraph.add((rdfNode, definition_node, rdflib.Literal(
            "A subsumption category for glycans described by the arrangement of monosaccharides and the glycosidic bonds linking them, and with partial or complete linkage position or anomeric configuration indicated.")))
        outputGraph.add((rdfNode, rdflib.namespace.RDFS.seeAlso,
                         rdflib.URIRef("http://www.glycoinfo.org/glyco/owl/relation#Linkage_defined_saccharide")))

        # Glycan class under OWL things
        outputGraph.add((ns["GNO_00000001"], rdflib.namespace.RDF.type, rdflib.OWL.Class))
        outputGraph.add((ns["GNO_00000001"], rdflib.namespace.RDFS.label, rdflib.Literal("glycan")))
        outputGraph.add((ns["GNO_00000001"], definition_node,
                         rdflib.Literal("A compound consisting of monosaccharides linked by glycosidic bonds.")))
        outputGraph.add((ns["GNO_00000001"], rdflib.namespace.RDFS.seeAlso,
                         rdflib.URIRef("http://purl.obolibrary.org/obo/CHEBI_50699")))
        outputGraph.add((ns["GNO_00000001"], rdflib.namespace.RDFS.seeAlso,
                         rdflib.URIRef("http://purl.obolibrary.org/obo/CHEBI_18154")))

        for n in self._nodes.values():
            if n._nodeType != "molecularweight":
                rdfNode = ns["GNO_%s" % n.getID()]
            else:
                try:
                    id = self.massiddict[float(n.getID())]
                except KeyError:
                    mass = n.getID()
                    id = str(int(max(self.massiddict.values())) + 1)
                    self.massiddict[float(mass)] = id
                    self.newMass = True
                rdfNode = ns["GNO_%s" % id]
            outputGraph.add((rdfNode, rdflib.namespace.RDF.type, rdflib.OWL.Class))

            if n._nodeType:
                # outputGraph.add((rdfNode, rdflib.namespace.RDFS.label, ns2[n._nodeType]))
                outputGraph.add((rdfNode, has_subsumption_level_node, subsumptionLevel[n._nodeType.lower()]))
            else:
                raise ValueError

            if n._nodeType != "molecularweight":
                outputGraph.add((rdfNode, has_glytoucan_id_node, rdflib.Literal(n.getID())))
                outputGraph.add((rdfNode, has_glytoucan_link_node, nsgtc[n.getID()]))
            else:
                outputGraph.add((rdfNode, rdflib.namespace.RDFS.subClassOf, ns["GNO_00000001"]))
                outputGraph.add((rdfNode, rdflib.namespace.RDFS.label,
                                 rdflib.Literal("glycan of molecular weight %s Da." % n.getID())))
                outputGraph.add((rdfNode, definition_node,
                                 rdflib.Literal(
                                     "A glycan characterized by underivitized molecular weight of %s Daltons" % n.getID())))

        for l in self.allRelationship():
            if l.getEdgeType() == "subsumes":
                if l._sideA._nodeType == "molecularweight":
                    id = self.massiddict[float(l._sideA.getID())]
                    n1 = ns["GNO_%s" % id]
                else:
                    n1 = ns["GNO_%s" % l._sideA.getID()]
                n2 = ns["GNO_%s" % l._sideB.getID()]
                outputGraph.add((n2, rdflib.namespace.RDFS.subClassOf, n1))

        if self.newMass:
            self.overwritemasslookuptable()

        return outputGraph

    def OWLoutput(self):
        return self.OWLgenerate().serialize(format="pretty-xml")


class dumpFile:
    nodes = {}
    edges = {}

    def readfile(self, dumpfilepath):
        f = open(dumpfilepath)
        raw_data = {}
        for l in f:
            l = l.strip()
            if l.startswith("#"):
                lineInfo = "node"  # node type none
                if l.startswith("# NODES"):
                    lineInfo = "node"
                    mass = float(re.compile("\d{2,5}\.\d{1,2}").findall(l)[0])
                    mass = "%.2f" % mass
                    content = {"nodes": {}, "edges": {}}
                    raw_data[mass] = content
                elif l.startswith("# EDGES"):
                    lineInfo = "edge"
                else:
                    lineInfo = "none"
            else:
                if l.startswith("WARNING"):
                    continue

                if lineInfo == "none":
                    continue
                elif lineInfo == "node":
                    nodeacc = l.split()[0]
                    nodetype = l.split()[2]
                    nodetype = nodetype.rstrip("*")
                    content["nodes"][nodeacc] = nodetype
                elif lineInfo == "edge":
                    to = re.compile("G\d{5}\w{2}").findall(l)
                    fromx = to.pop(0)
                    if to:
                        content["edges"][fromx] = to
                else:
                    raise RuntimeError
        return raw_data


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

        g = SubsumptionGraph(*sys.argv[1:], verbose=verbose)

    elif cmd == "writeowl":
        def write_owl():
            r = relationship()
            df = dumpFile()
            # "gnome_subsumption_raw.txt"
            ifn = sys.argv[1] # path to input file
            data = df.readfile(ifn)

            for mass, v in data.items():
                mass = str(float(mass))
                # print mass
                nodes = v["nodes"]
                edges = v["edges"]

                r.addNode(mass, nodetype="molecularweight")
                notParents = []
                for son in edges.values():
                    notParents += son

                for n, t in nodes.items():
                    r.addNode(n, nodetype=t)

                for p, cs in edges.items():
                    # print mass
                    for c in cs:
                        # print p,c
                        r.connect(r.getNode(p), r.getNode(c), "subsumes")

                for potentialr in (set(nodes.keys()) - set(notParents)):
                    r.connect(r.getNode(mass), r.getNode(potentialr), "subsumes")
                # print nodes
                # print edges
            if len(sys.argv) > 2:
                ofn = sys.argv[2]
            else:
                ofn = "GNOme_temp.owl"
            f = open(ofn, "w")
            s = r.OWLoutput()
            f.write(s)
        write_owl()

    else:

        print >> sys.stderr, "Bad command: %s" % (cmd,)
        sys.exit(1)
