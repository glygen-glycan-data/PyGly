import copy
import re
import sys
import os, os.path, urllib
from collections import defaultdict

import rdflib


class GNOme(object):
    version = "1.1.2"
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
        self.ns['gno'] = rdflib.Namespace('http://purl.obolibrary.org/obo/GNO_')
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
                attr[plab] = olab
            else:
                olab = self.accession(o)
                if plab not in attr:
		                attr[plab] = []
                attr[plab].append(olab)
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
        writer = OWLWriter()
        writer.write(handle, self.gnome)

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
		    # cluster[acc]['level'] = 'Saccharide*'
                continue
            if self.any_parent_pos(gly):
                if not level:
                    cluster[acc]['level'] = 'Saccharide*'
                elif level != "Saccharide":
                    self.warning(
                        "annotation inferred level %s for %s != computed level Saccharide (parent_pos)" % (level, acc),
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
            if self.monosaccharide_count(gly) == 1 and self.any_ring(gly):
                if not level:
                    cluster[acc]['level'] = 'Topology*'
                elif level != "Topology":
                    self.warning("annotation inferred level %s for %s != computed level Topology" % (level, acc), 1)
                    # cluster[acc]['level'] = 'Topology*'
                continue
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
            if l.parent_pos() != None and l.parent_pos() != set([l.parent().ring_start()]):
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


    # Dump file parsing part starts here
    raw_data = {}
    dumpfilepath = ""
    warnings = defaultdict(set)
    warningsbytype = None

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
                    content = {"nodes": {}, "edges": {}}
                    raw_data[mass] = content
                elif l.startswith("# EDGES"):
                    lineInfo = "edge"
                else:
                    lineInfo = "none"
            else:

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
        self.raw_data = raw_data

        self.allnodestype = {}
        self.alledges = {}
        for component in raw_data.values():
            self.allnodestype.update(copy.deepcopy(component["nodes"]))
            self.alledges.update(copy.deepcopy(component["edges"]))

        allmass = list()
        for m in self.raw_data.keys():
            self.allnodestype[m] = "molecular weight"
            allmass.append(m)
            top = set(raw_data[m]["nodes"].keys())
            for chilren in raw_data[m]["edges"].values():
                top = top - set(chilren)
            top = list(top)
            self.alledges[m] = top

        self.allnodestype["00000001"] = "glycan"
        self.alledges["00000001"] = allmass

        return raw_data

    def root(self):
        return "00000001"

    def nodes(self):
        for n in self.allnodestype.keys():
            yield n

    def edges(self):
        for n in self.nodes():
            if n in self.alledges:
                for c in self.alledges[n]:
                    yield (n, c)

    def parents(self, accession):
        for p, c in self.alledges.items():
            if accession in c:
                yield p

    def ancestors(self, accession):
        anc = set()
        for p in self.parents(accession):
            anc.add(p)
            anc.update(self.ancestors(p))
        return anc

    def children(self, accession):
        for c in self.alledges.get(accession, []):
            yield c

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
        return self.allnodestype[accession]

    def islevel(self, accession, level):
        return self.level(accession) == level

    def ismolecularweight(self, accession):
        return self.islevel(accession, 'molecular weight')

    def isbasecomposition(self, accession):
        return self.islevel(accession, 'BaseComposition')

    def iscomposition(self, accession):
        return self.islevel(accession, 'Composition')

    def istopology(self, accession):
        return self.islevel(accession, 'Topology')

    def issaccharide(self, accession):
        return self.islevel(accession, 'Saccharide')

    def get_basecomposition(self, accession):
        for n, t in self.allnodestype.items():
            if self.isbasecomposition(n) and accession in self.descendants(n):
                yield n

    def get_composition(self, accession):
        for n, t in self.allnodestype.items():
            if self.iscomposition(n) and accession in self.descendants(n):
                yield n

    def get_topology(self, accession):
        for n, t in self.allnodestype.items():
            if self.istopology(n) and accession in self.descendants(n):
                yield n

    def has_basecomposition(self, accession):
        assert self.isbasecomposition(accession)
        return self.descendants(accession)

    def has_composition(self, accession):
        assert self.iscomposition(accession)
        return self.descendants(accession)

    def has_topology(self, accession):
        assert self.istopology(accession)
        return self.descendants(accession)

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


from rdflib import URIRef, Namespace


class OWLWriter():
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

    gno = "http://purl.obolibrary.org/obo/"
    gtcs = "http://glytoucan.org/Structures/Glycans/"
    iao = "http://purl.obolibrary.org/obo/"
    gtco = "http://www.glytoucan.org/glyco/owl/glytoucan#"
    rocs = "http://www.glycoinfo.org/glyco/owl/relation#"

    glycan_class = 1
    definition_class = "IAO_0000115"
    subsumption_level_class = 11
    subsumption_level_annotation_property = 21
    glytoucan_id_annotation_property = 22
    glytoucan_link_annotation_property = 23

    subsumption_level = {

        "molecularweight": {"id": 12,
                            "label": "subsumption category molecular weight",
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
                            "label": "subsumption category basecomposition",
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
                        "label": "subsumption category composition",
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
                     "label": "subsumption category topology",
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
                       "label": "subsumption category saccharide",
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
        dc = rdflib.namespace.DC

        Literal = rdflib.Literal

        gno = Namespace(self.gno)
        gtcs = Namespace(self.gtcs)
        iao = Namespace(self.iao)
        gtco = Namespace(self.gtco)
        rocs = Namespace(self.rocs)

        outputGraph.bind("owl", owl)
        outputGraph.bind("gno", gno)
        outputGraph.bind("obo", iao)
        outputGraph.bind("dc", dc)
        outputGraph.bind("rocs", rocs)

        root = rdflib.URIRef("http://purl.obolibrary.org/obo/gno.owl")

        # Add ontology
        outputGraph.add((root, rdf.type, owl.Ontology))

        # Copyright
        outputGraph.add(
            (root, dc.license, URIRef("http://creativecommons.org/licenses/by/4.0/")))
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
            if n._nodeType != "molecularweight":
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

            if n._nodeType != "molecularweight":
                outputGraph.add((rdfNode, has_glytoucan_id_node, Literal(n.getID())))
                outputGraph.add((rdfNode, has_glytoucan_link_node, gtcs[n.getID()]))
            else:
                outputGraph.add((rdfNode, rdfs.subClassOf, self.gnouri(self.glycan_class)))
                outputGraph.add((rdfNode, rdfs.label,
                                 Literal("glycan of molecular weight %s Da." % n.getID())))
                outputGraph.add((rdfNode, definition,
                                 Literal(
                                     "A glycan characterized by underivitized molecular weight of %s Daltons" % n.getID())))

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
        r = OWLWriter()

        # "../smw/glycandata/data/gnome_subsumption_raw.txt"
        ifn = sys.argv[1]  # path to input file
        df = dumpFile(ifn)

        for mass in df.getmass():
            mass = "%.2f" % float(mass)
            # print mass
            nodes = df.getnodesbymass(mass)

            r.addNode(mass, nodetype="molecularweight")

            for n in nodes:
                r.addNode(n, nodetype=df.getnodetype(n))

            for n in nodes:
                if df.getchild(n):
                    for c in df.getchild(n):
                        r.connect(r.getNode(n), r.getNode(c), "subsumes")

            for n in df.getnodessubsumedbymass(mass):
                r.connect(r.getNode(mass), r.getNode(n), "subsumes")

        if len(sys.argv) > 2:
            ofn = sys.argv[2]
        else:
            ofn = "GNOme_temp.owl"
        f = open(ofn, "w")
        s = r.write(f, r.make_graph())
        f.close()


    else:

        print >> sys.stderr, "Bad command: %s" % (cmd,)
        sys.exit(1)
