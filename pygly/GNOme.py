import copy
import re
import sys
import ssl
import os, os.path, urllib, time
import psutil
from collections import defaultdict
import datetime

import rdflib
from rdflib.namespace import XSD
import json

def printmem():
    process = psutil.Process()
    mi = process.memory_info()
    res = float(mi.rss)
    vir = float(mi.vms)
    print("# Memory: %.2fGB (virtual) %.2fGB (resident)"%(res/1024**3,vir/1024**3))

class GNOmeAPI(object):

    # Base class for GNOme and subsumption API supported by both the
    # released OWL file and the so-called Subsumption Graph raw dump
    # file.

    # Subsumption types - Ugh!
    class EdgeType:
        Any = None
        Strict = 1
        Other = 2

    # All nodes
    def nodes(self):
        raise NotImplemented("GNOme API method nodes not implemented")

    # Remove a node
    def delete_node(self, accession):
        raise NotImplemented("GNOme API method delete_node not implemented")

    # Child nodes of a node
    def children(self, accession, edgetype=EdgeType.Strict):
        raise NotImplemented("GNOme API method children not implemented")

    # Parent nodes of a node
    def parents(self, accession, edgetype=EdgeType.Strict):
        raise NotImplemented("GNOme API method parents not implemented")

    # Set the parents of a node
    def set_parents(self, accession, parents, edgetype=EdgeType.Strict):
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
    def edges(self,edgetype=EdgeType.Strict):
        for n in self.nodes():
            for c in self.children(n, edgetype):
                yield n,c

    # Descendants of a node
    def descendants(self, accession, edgetype=EdgeType.Strict):
        desc = set()
        for c in self.children(accession, edgetype):
            desc.add(c)
            desc.update(self.descendants(c, edgetype))
        return desc

    # Ancestors of a node
    def ancestors(self, accession, edgetype=EdgeType.Strict):
        anc = set()
        for p in self.parents(accession, edgetype):
            anc.add(p)
            anc.update(self.ancestors(p, edgetype))
        return anc

    # Whether a node is a leaf
    def isleaf(self, accession, edgetype=EdgeType.Strict):
        for ch in self.children(accession, edgetype):
            return False
        return True

    # Whether a node is a root node
    def isroot(self,accession):
        for pt in self.parents(accession, edgetype=EdgeType.Strict):
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

    def restrict(self, restriction, **kw):

        # Update the restriction set to include "landmark" topology, composition, etc nodes
        keep = set(restriction)
        keep.add(self.root())
        for acc in restriction:
            keep.add(self.get_archetype(acc))
            keep.add(self.get_topology(acc))
            keep.add(self.get_composition(acc))
            keep.add(self.get_basecomposition(acc))
            keep.add(self.get_molecularweight(acc))
            if kw.get('allstructanc',False) or kw.get('alltopoanc',False):
                for anc in self.ancestors(acc):
                    if (kw.get('allstructanc',False) and self.issaccharide(anc)) or \
                       (kw.get('alltopoanc',False) and self.istopology(anc)):
                        keep.add(anc)
        if None in keep:
            keep.remove(None)

        for acc in keep:
            if self.get_archetype(acc) not in keep:
                self.delete_predicate(acc,"00000036")

        for acc in restriction:
             self.gnome.add((self.uri("gno:" + acc), self.uri("gno:00000025"), rdflib.Literal("true",datatype=rdflib.XSD.boolean)))

        # find all ancestors of each kept node
        parents = defaultdict(set)
        for acc in keep:
            for anc in self.ancestors(acc):
                if anc in keep:
                    parents[acc].add(anc)

        # find all ancestors of each kept node
        allparents = defaultdict(set)
        for acc in keep:
            for anc in self.ancestors(acc,GNOmeAPI.EdgeType.Any):
                if anc in keep:
                    allparents[acc].add(anc)

        # then eliminate shortcuts
        toremove = set()
        for n1 in keep:
            for n2 in parents[n1]:
                for n3 in parents[n2]:
                    if n3 in parents[n1]:
                        toremove.add((n1, n3))
        for n1, n3 in toremove:
            parents[n1].remove(n3)

        # then eliminate shortcuts
        toremove = set()
        for n1 in keep:
            for n2 in allparents[n1]:
                for n3 in allparents[n2]:
                    if n3 in allparents[n1]:
                        toremove.add((n1, n3))
        for n1, n3 in toremove:
            allparents[n1].remove(n3)

        for n in self.nodes():
            if n not in keep:
                self.delete_node(n)
            else:
                self.set_parents(n,parents[n],GNOmeAPI.EdgeType.Strict)
                self.set_parents(n,allparents[n]-parents[n],GNOmeAPI.EdgeType.Other)

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
        self.version = None
        self.versiondate = None
        # for s,p,o in self.triples(None,None,None):
        #     print(s,p,o)
        for s,p,o in self.triples(None,"owl:versionIRI",None):
            versionurl = str(o)
            break
        if versionurl:
            self.versiondate = versionurl.split('/')[-2]
        for s,p,o in self.triples(None,"owl:versionInfo",None):
            self.version = str(o).lower()
            break

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
                return str(o)
        for ns in self.ns.values():
            if ns.startswith('http://') and uri.startswith(ns):
                return uri[len(ns):]
        return str(uri)

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

    def parents(self, accession, edgetype=GNOmeAPI.EdgeType.Strict):
        uri = "gno:%s" % (accession,)
        if edgetype == GNOmeAPI.EdgeType.Strict:
            for s, p, o in self.triples(uri, 'rdfs:subClassOf', None):
                yield self.accession(o)
        elif edgetype == GNOmeAPI.EdgeType.Any:
            for s, p, o in self.triples(uri, 'gno:00000024', None):
                yield self.accession(o)
        elif edgetype == GNOmeAPI.EdgeType.Other:
            items = set()
            for s, p, o in self.triples(uri, 'gno:00000024', None):
                items.add(self.accession(o))
            for s, p, o in self.triples(uri, 'rdfs:subClassOf', None):
                items.remove(self.accession(o))
            for acc in items:
                yield acc

    def children(self, accession, edgetype=GNOmeAPI.EdgeType.Strict):
        uri = "gno:%s" % (accession,)
        if edgetype == GNOmeAPI.EdgeType.Strict:
            for s, p, o in self.triples(None, 'rdfs:subClassOf', uri):
                yield self.accession(s)
        elif edgetype == GNOmeAPI.EdgeType.Any:
            for s, p, o in self.triples(None, 'gno:00000024', uri):
                yield self.accession(s)
        elif edgetype == GNOmeAPI.EdgeType.Other:
            items = set()
            for s, p, o in self.triples(None, 'gno:00000024', uri):
                items.add(self.accession(s))
            for s, p, o in self.triples(None, 'rdfs:subClassOf', uri):
                items.remove(self.accession(s))
            for acc in items:
                yield acc

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

    def get_archetype(self, accession):
        uri = "gno:%s" % (accession,)
        for s, p, o in self.triples(uri, "gno:00000036", None):
            return self.label(o)

    def get_restrictionmember(self, accession):
        uri = "gno:%s" % (accession,)
        for s, p, o in self.triples(uri, "gno:00000025", None):
            return o.eq(True)
        return False

    def isarchetype(self, accession):
        return self.get_archetype(accession) == accession

    def delete_node(self, accession):
        self.gnome.remove((self.uri("gno:" + accession), None, None))

    def delete_predicate(self, accession, predicate):
        self.gnome.remove((self.uri("gno:" + accession), self.uri("gno:" + predicate), None))

    def set_parents(self, accession, parents, edgetype=GNOmeAPI.EdgeType.Strict):
        assert edgetype in (GNOmeAPI.EdgeType.Strict,GNOmeAPI.EdgeType.Other)
        strict = set(self.parents(accession,GNOmeAPI.EdgeType.Strict))
        other = set(self.parents(accession,GNOmeAPI.EdgeType.Other))
        self.gnome.remove((self.uri("gno:" + accession), self.uri("rdfs:subClassOf"), None))
        self.gnome.remove((self.uri("gno:" + accession), self.uri("gno:00000024"), None))
        if edgetype == GNOmeAPI.EdgeType.Strict:
            strict = set(parents)
        else:
            other = set(parents)
        for parent in strict:
            self.gnome.add((self.uri("gno:" + accession), self.uri("rdfs:subClassOf"), self.uri("gno:" + parent)))
            self.gnome.add((self.uri("gno:" + accession), self.uri("gno:00000024"), self.uri("gno:" + parent)))
        for parent in other:
            self.gnome.add((self.uri("gno:" + accession), self.uri("gno:00000024"), self.uri("gno:" + parent)))

    def write(self, handle):
        writer = OWLWriter()
        writer.write(handle, self.gnome)

    def dump(self):
        for acc in sorted(g.nodes()):
            print(acc)
            for k, v in sorted(g.attributes(acc).items()):
                print("  %s: %s" % (k, v))

    def get_cb_button_str(self, acc):
        for s, p, o in self.triples(None, "gno:00000101", None):
            if acc == self.accession(s):
                return o

    # ics2dp: IUPAC composition string to dictionary pattern(regex)
    ics2dp = re.compile(r"([\D\+]{1,13})(\d{1,3})")
    def iupac_composition_str_to_dict(self, s):
        res = {}
        for p in self.ics2dp.findall(s):
            res[p[0]] = int(p[1])

        additional_xxx = 0
        additional_xxx += res.get('HexNAc+aldi',0)
        additional_xxx += res.get('Hex+aldi',0)
        additional_xxx += res.get('dHex+aldi',0)
 
        for mmm, count in list(res.items()):
            if "+aldi" in mmm:
                del res[mmm]

        if additional_xxx > 0:
            if "Xxx" not in res:
                res["Xxx"] = additional_xxx
            else:
                res["Xxx"] = res["Xxx"]+additional_xxx
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

    def get_structure_characterization_score(self, acc):
        for s, p, o in self.triples(None, "gno:00000102", None):
            if acc == self.accession(s):
                return o

    def all_structure_characterization_score(self):
        res = {}
        for s, p, o in self.triples(None, "gno:00000102", None):
            acc = self.accession(s)
            res[acc] = o
        return res

    def toViewerData(self, output_file_path):
        assert(self.version)
        res = {'__VERSION__': self.version}
        cbbutton = self.all_cbbutton()
        # byonic = self.all_Byonic()
        syms = self.all_synonym()
        scores = self.all_structure_characterization_score()

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

            arch = self.get_archetype(n)
            isarch = self.isarchetype(n)

            if n not in cbbutton:
                continue

            children = sorted(self.children(n,GNOmeAPI.EdgeType.Any))
            # parents = sorted(self.parents(n,GNOmeAPI.EdgeType.Any))
            otherchildren = sorted(self.children(n,GNOmeAPI.EdgeType.Other))
            edgetype = dict()
            for acc in otherchildren:
                edgetype[acc] = 'other'

            res[n] = {
                "level": t, "children": children, "edgetype": edgetype, "count": cbbutton[n], "score": int(scores[n]), "archetype": arch
            }
            if len(children) == 0:
                del res[n]['children']
            if len(otherchildren) == 0:
                del res[n]['edgetype']
            if arch == None:
                del res[n]['archetype']

            if n in syms:
                res[n]['syms'] = sorted(syms[n])

            if self.get_restrictionmember(n):
                res[n]['inrestriction'] = True

        json.dump(res, open(output_file_path, 'w'), indent=2, separators=(',',': '), sort_keys=True)
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




from . alignment import GlycanSubsumption, GlycanEqual, GlycanEqualWithWURCSCheck
from . Monosaccharide import Anomer, Substituent
from . GlycanResource import GlyTouCan, GlyCosmos
from . manipulation import Topology, Composition, BaseComposition, Archetype, WURCSArchetype, WURCSManipulation, RemoveAlditol


class SubsumptionGraph(GNOmeAPI):
    def __init__(self, *args, **kwargs):
        pass

    def compute(self, *args, **kwargs):
        self.gtc = GlyTouCan(usecache=kwargs.get('usecache',False))
        self.gco = GlyCosmos(usecache=kwargs.get('usecache',False))
        self.subsumption = GlycanSubsumption()
        self.geq = GlycanEqual()
        self.geqwwc = GlycanEqualWithWURCSCheck()
        self.archetype = Archetype()
        self.wurcs_archetype = WURCSArchetype()
        self.remove_aldi = RemoveAlditol()
        self.topology = Topology()
        self.composition = Composition()
        self.basecomposition = BaseComposition()
        self.verbose = kwargs.get('verbose', 0)
        self.score = IncompleteScore()

        # We establish the following sets
        # 1. GlyTouCan accessions (self.allacc)
        # 2. Archived GlyTouCan accession (some with replacements) (self.replace)
        # 3. GlyCosmos (validated) accessions (self.allgco)

        self.allacc = set(self.gtc.allaccessions())
        self.allgco = set(self.gco.allaccessions())
        self.replace = self.gco.replace()
        for acc in self.replace:
            try:
                self.allacc.remove(acc)
            except KeyError:
                pass

        printmem()

        masscluster = defaultdict(dict)
        extramasscluster = defaultdict(dict)
        clustermap = defaultdict(set)
        extraclustermap = defaultdict(set)
        if len(args) > 0:
            argmass = []
            for a in args:
                if '-' in a:
                    low,high = [ s.strip() for s in a.split('-') ]
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
            for glyacc in sorted(self.allacc):
                mass = self.gco.getmass(glyacc) 
                if not mass:
                    mass = self.gtc.getmass(glyacc)
                mass1 = self.gco.umw(accession=glyacc,format='wurcs')
                if not mass1:
                    mass1 = self.gtc.umw(accession=glyacc,format='wurcs')
                if not mass and not mass1:
                    self.warning("mass could not be determined for %s" % (glyacc), 5)
                    continue
                elif mass and mass1 and abs(mass - mass1) > 0.0001:
                    self.warning("mass inconsistency for %s: %f vs %f" % (glyacc, mass, mass1), 2)
                if mass1:
                    mass = mass1
                if not mass:
                    continue
                rmass = str(round(mass, 2))
                gly = self.gco.getGlycan(glyacc,format='wurcs')
                if not gly:
                    gly = self.gtc.getGlycan(glyacc,format='wurcs')
                if not gly:
                    continue
                if glyacc in clustermap:
                    for rmi in clustermap[glyacc]:
                        if glyacc in masscluster[rmi]:
                            del masscluster[rmi][glyacc]
                    clustermap[glyacc].add(rmass)
                    continue
                for low,high in argmass:
                    if float(low) <= float(rmass) <= float(high):
                        masscluster[rmass][glyacc] = dict(accession=glyacc,glycan=gly)
                        clustermap[glyacc].add(rmass)
                        break
                if gly and gly.has_alditol_root():
                    gly1 = self.remove_aldi(gly)
                    rmass_minus_aldi = str(round(mass-2.015650070,2))
                    if glyacc in extraclustermap:
                        for rmi in extraclustermap[glyacc]:
                            if glyacc in extramasscluster[rmi]:
                                del extramasscluster[rmi][glyacc]
                        extraclustermap[glyacc].add(rmass_minus_aldi)
                        continue
                    for low,high in argmass:
                        if float(low) <= float(rmass_minus_aldi) <= float(high):
                            extramasscluster[rmass_minus_aldi][glyacc] = dict(accession=glyacc,mass=mass,rmass=rmass,mass1=mass-2.015650070,rmass1=rmass_minus_aldi,delta=+2.015650070,glycan=gly1)
                            extraclustermap[glyacc].add(rmass_minus_aldi)
        else:
            for glyacc in sorted(self.allacc):
                mass = self.gco.getmass(glyacc)
                if not mass:
                    mass = self.gtc.getmass(glyacc)
                mass1 = self.gco.umw(accession=glyacc,format='wurcs')
                if not mass1:
                    mass1 = self.gtc.umw(accession=glyacc,format='wurcs')
                if not mass and not mass1:
                    self.warning("mass could not be determined for %s" % (glyacc), 5)
                    continue
                elif mass and mass1 and abs(mass - mass1) > 0.0001:
                    self.warning("mass inconsistency for %s: %f vs %f" % (glyacc, mass, mass1), 2)
                if mass1:
                    mass = mass1
                if not mass:
                    continue
                rmass = str(round(mass, 2))
                gly = self.gco.getGlycan(glyacc,format='wurcs')
                if not gly:
                    gly = self.gtc.getGlycan(glyacc,format='wurcs')
                if glyacc in clustermap:
                    for rmi in clustermap[glyacc]:
                        if glyacc in masscluster[rmi]:
                            del masscluster[rmi][glyacc]
                    clustermap[glyacc].add(rmass)
                    continue
                masscluster[rmass][glyacc] = dict(accession=glyacc,glycan=gly)
                clustermap[glyacc].add(rmass)
                if gly and gly.has_alditol_root():
                    gly1 = self.remove_aldi(gly)
                    rmass_minus_aldi = str(round(mass-2.015650070,2))
                    if glyacc in extraclustermap:
                        for rmi in extraclustermap[glyacc]:
                            if glyacc in extramasscluster[rmi]:
                                del extramasscluster[rmi][glyacc]
                        extraclustermap[glyacc].add(rmass_minus_aldi)
                        continue
                    extramasscluster[rmass_minus_aldi][glyacc] = dict(accession=glyacc,mass=mass,rmass=rmass,mass1=mass-2.015650070,rmass1=rmass_minus_aldi,delta=+2.015650070,glycan=gly1)
                    extraclustermap[glyacc].add(rmass_minus_aldi)

        for rmass,cluster in masscluster.items():
            for acc in list(cluster):
                topo = self.gtc.gettopo(acc)
                if topo in self.replace:
                    topo = self.replace[topo]
                if topo and topo not in cluster:
                    masscluster[rmass][topo] = dict(accession=topo)
                    for rmi in clustermap.get(topo,[]):
                        del masscluster[rmi][topo]
                comp = self.gtc.getcomp(acc)
                if comp in self.replace:
                    comp = self.replace[comp]
                if comp and comp not in cluster:
                    masscluster[rmass][comp] = dict(accession=comp)
                    for rmi in clustermap.get(comp,[]):
                        del masscluster[rmi][comp]
                bcomp = self.gtc.getbasecomp(acc)
                if bcomp in self.replace:
                    bcomp = self.replace[bcomp]
                if bcomp and bcomp not in cluster:
                    masscluster[rmass][bcomp] = dict(accession=bcomp)
                    for rmi in clustermap.get(bcomp,[]):
                        del masscluster[rmi][bcomp]

        printmem()

        for rmass, cluster in sorted(masscluster.items(),key=lambda t: float(t[0])):
            self.compute_component(rmass, masscluster, extramasscluster)
            printmem()

    def warning(self, msg, level):
        if self.verbose >= level:
            print("# WARNING:%d - %s" % (level, msg))
            sys.stdout.flush()

    def compute_component(self, rmass, masscluster, extramasscluster):
        cluster = masscluster[rmass];
        extracluster = extramasscluster[rmass]
        combinedcluster = dict()
        combinedcluster.update(cluster.items())
        combinedcluster.update(extracluster.items())
        start = time.time()

        print("# START %s - %d accessions in molecular weight cluster for %s" % (time.ctime(), len(cluster), rmass))
        sys.stdout.flush()

        badparse = 0
        total = len(cluster)
        allgly = dict()
        for acc in sorted(cluster):
            if acc in self.replace:
                repl = self.replace[acc]
                if repl == None:
                    self.warning("Glycan accession archived: %s"%(acc,),2)
                else:
                    self.warning("Glycan accession %s archived and replaced by %s"%(acc,repl,),2)
                continue
            gly = cluster[acc].get('glycan')
            if not gly:
                gly = self.gco.getGlycan(acc,format='wurcs')
            if not gly:
                gly = self.gtc.getGlycan(acc,format='wurcs')
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
            if acc not in self.allgco:
                self.warning("Glycan not validated by GlyCosmos: %s"%(acc,),3)
                # continue
            gcomass = self.gco.getmass(acc)
            gtcmass = self.gtc.getmass(acc)
            mygcomass = self.gco.umw(accession=acc,format='wurcs')
            mygtcmass = self.gtc.umw(accession=acc,format='wurcs')
            if gcomass and abs(float(rmass)-gcomass) <= 0.01:
                cluster[acc]['mass'] = gcomass
            elif gtcmass and abs(float(rmass)-gtcmass) <= 0.01:
                cluster[acc]['mass'] = gtcmass
            elif mygcomass and abs(float(rmass)-mygcomass) <= 0.01:
                cluster[acc]['mass'] = mygcomass
            elif mygtcmass and abs(float(rmass)-mygtcmass) <= 0.01:
                cluster[acc]['mass'] = mygtcmass
            else:
                cluster[acc]['mass'] = float(rmass)

        clusteracc = set(map(lambda t: t[0], filter(lambda t: t[1].get('glycan'), cluster.items())))
        combclusteracc = set(map(lambda t: t[0], filter(lambda t: t[1].get('glycan'), combinedcluster.items())))

        outedges = defaultdict(set)
        inedges = defaultdict(set)

        self.warning("Computation of subsumption relationships started",5)
        for acc1 in sorted(combclusteracc):
            gly1 = combinedcluster[acc1]['glycan']
            for acc2 in sorted(combclusteracc):
                gly2 = combinedcluster[acc2]['glycan']
                if acc1 != acc2 and not (acc1 in cluster and acc2 in extracluster):
                    self.warning("%s <?= %s"%(acc1,acc2),5)
                    if self.subsumption.leq(gly1, gly2):
                        iseq = self.geq.eq(gly1, gly2)
                        if not iseq or acc2 < acc1 or (acc1 in extracluster and acc2 in cluster):
                            outedges[acc2].add(acc1)
                            inedges[acc1].add(acc2)
                        if iseq and (acc1 in cluster and acc2 in cluster) and self.geqwwc.eq(gly1, gly2) and acc1 < acc2:
                            self.warning("Potential WURCS canonicalization issue: %s == %s"%(acc1,acc2),1)
        self.warning("Computation of subsumption relationships done",5)

        # Check GlyTouCan topology, composition, basecomposition, and
        # molecular weight annotations with respect to computed
        # subsumption relationships and computed molecular weight

        for acc in sorted(clusteracc):

            topo = self.gtc.gettopo(acc)
            if topo in self.replace:
                topo = self.replace[topo]
            if topo:
                if topo not in cluster:
                    self.warning("annotated topology %s of %s is not in %s rounded mass cluster" % (topo, acc, rmass), 1)
                elif not cluster[topo].get('glycan'):
                    self.warning("annotated topology %s of %s cannot be parsed" % (topo, acc), 1)
                elif acc not in outedges[topo] and acc != topo:
                    self.warning("annotated topology %s does not subsume %s" % (topo, acc), 1)

            comp = self.gtc.getcomp(acc)
            if comp in self.replace:
                comp = self.replace[comp]
            if comp:
                if comp not in cluster:
                    self.warning("annotated composition %s of %s is not in %s rounded mass cluster" % (comp, acc, rmass), 1)
                elif not cluster[comp].get('glycan'):
                    self.warning("annotated composition %s of %s cannot be parsed" % (comp, acc), 1)
                elif acc not in outedges[comp] and acc != comp:
                    self.warning("annotated composition %s does not subsume %s" % (comp, acc), 1)

            bcomp = self.gtc.getbasecomp(acc)
            if bcomp in self.replace:
                bcomp = self.replace[bcomp]
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
            if bcomp in self.replace:
                bcomp = self.replace[bcomp]
            comp = self.gtc.getcomp(acc)
            if comp in self.replace:
                comp = self.replace[comp]
            topo = self.gtc.gettopo(acc)
            if topo in self.replace:
                topo = self.replace[topo]

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

                # if False and level1 == "Topology" and acc == "G65022XW":
                #     print self.geq.eq(gly1,self.topology(gly))
                #     print "topology(%s)"%(acc,)
                #     print self.topology(gly).glycoct()
                #     print acc1
                #     print gly1.glycoct()

                # if False and level1 == "Composition" and acc == "G76453ML":
                #     print self.geq.eq(gly1,self.composition(gly))
                #     print "composition(%s)"%(acc,)
                #     print self.composition(gly).glycoct()
                #     print acc1
                #     print gly1.glycoct()
                
                self.warning("Checking topo, comp, bcomp relationships for %s: %s (%s)"%(acc,acc1,level1),5)
                if level1 == "BaseComposition" and self.geq.eq(gly1,self.basecomposition(gly)) and acc1 not in self.replace:
                    bcomp.add(acc1)
                elif level1 == "Composition" and self.geq.eq(gly1,self.composition(gly)) and acc1 not in self.replace:
                    comp.add(acc1)
                elif level1 == "Topology" and self.geq.eq(gly1,self.topology(gly)) and acc1 not in self.replace:
                    topo.add(acc1)

            if len(topo) > 1:
                self.warning("multiple topologies %s for %s"%(", ".join(topo), acc), 1)
                if g.get('topo') in topo:
                    # Take GTC one if present
                    topo = g.get('topo')
                else:
                    # Take the one subsumed by all the others...
                    topo1=None
                    for acc in topo-set(self.replace):
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
                    for acc in comp-set(self.replace):
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
                    for acc in bcomp-set(self.replace):
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

        for acc0 in clusteracc:
            gly0 = cluster[acc0]['glycan']
            # if acc0 == "G00026MO":
            #     print(acc0)
            #     print(gly0.glycoct())
            for acc1 in clusteracc:
                gly1 = cluster[acc1]['glycan']
                if gly1.has_root() and not gly1.repeated() and self.geq.eq(gly0,self.archetype(gly1)):
                    if self.wurcs_archetype.redend_mono(gly1.root().external_descriptor()) == gly0.root().external_descriptor() and (acc0 == acc1 or 'has_archetype' not in cluster[acc1]):
                        cluster[acc1]['has_archetype'] = acc0
            for acc1 in extracluster:
                gly1 = extracluster[acc1]['glycan']
                # if acc0 == "G00026MO" and acc1 == "G62167DK":
                #     print(acc1)
                #     print(gly1.glycoct())
                #     print(self.archetype(gly1).glycoct())
                if gly1.has_root() and not gly1.repeated() and self.geq.eq(gly0,self.archetype(gly1)):
                    if self.wurcs_archetype.redend_mono(gly1.root().external_descriptor()) == gly0.root().external_descriptor():
                        # extracluster[acc1]['has_archetype'] = acc0
                        if acc1 in masscluster[extracluster[acc1]['rmass']]:
                            masscluster[extracluster[acc1]['rmass']][acc1]['has_archetype'] = acc0 + "*"

        for acc in clusteracc:
            cluster[acc]['missing'] = None
            try:
                cluster[acc]['missing'] = self.score.score(cluster[acc]["glycan"])
            except LookupError:
                self.warning("Unable to compute missing rank for %s" % (acc), 1)

        if len(clusteracc) < 1:
            print("# DONE - Elapsed time %.2f sec." % (time.time() - start,))
            sys.stdout.flush()
            return

        print("# NODES - %d/%d + %d glycans in molecular weight cluster for %s" % (len(clusteracc), total, len(extracluster), rmass))
        for acc in sorted(clusteracc, key=lambda acc: (cluster[acc].get('level').rstrip('*'),acc) ):
            g = cluster[acc]
            print(acc,end=" ")
            print(g.get('mass'),end=" ")
            print(g.get('level'), g.get('topo'), g.get('comp'), g.get('bcomp'),end=" ")
            gly = g.get('glycan')
            extras = []
            if not gly.has_root():
                extras.append("COMP")
            if gly.undetermined() and gly.has_root():
                extras.append("UNDET")
            if gly.fully_determined():
                extras.append("FULL")
            if g.get('has_archetype'):
                extras.append("HASARCH:%s"%(g['has_archetype'],))
            for m,c in sorted(gly.iupac_composition(floating_substituents=True,aggregate_basecomposition=True).items()):
                if m != "Count" and c > 0:
                    extras.append("%s:%d"%(m,c))
            print(" ".join(extras))
        print("# ENDNODES - %d/%d + %d glycans in molecular weight cluster for %s" % (len(clusteracc), total, len(extracluster), rmass))
        sys.stdout.flush()

        prunedoutedges = self.prune_edges(outedges)
        prunedinedges = self.prune_edges(inedges)

        print("# TREE")
        for r in sorted(combclusteracc):
            if len(prunedinedges[r]) == 0:
                self.print_tree(combinedcluster, prunedoutedges, r)
        print("# ENDTREE")

        print("# EDGES")
        for n in sorted(clusteracc):
            print("%s:"%(n,),end=" ")
            print(" ".join(map(lambda s: s+"*" if "delta" in combinedcluster[s] else s,sorted(prunedoutedges[n]))))
        print("# ENDEDGES")
        sys.stdout.flush()

        print("# DONE - Elapsed time %.2f sec." % (time.time() - start,))
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
                l = l.parent().any_parent_link()
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
        
        print("%s%s%s" % (" " * indent, root, "*" if "delta" in cluster[root] else ""), cluster[root].get('level',"-"), cluster[root].get('missing','-'))
        for ch in sorted(edges[root]):
            self.print_tree(cluster, edges, ch, indent + 2)


    # Dump file parsing part starts here
    raw_data = {}
    dumpfilepath = ""
    warnings = defaultdict(set)
    warningsbytype = None
    monosaccharide_count = {}
    monosaccharide_count_pattern = re.compile(r"[\w\+]{1,13}:\d{1,3}")
    extra_properties = defaultdict(dict)

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
                    mass = float(re.compile(r"\d{2,6}\.\d{1,2}").findall(l)[0])
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
                    for cell in l.split()[6:]:
                        temp2 = self.monosaccharide_count_pattern.findall(cell)
                        if len(temp2) == 1:
                            temp2 = temp2[0].split(":")
                            iupac_comp_mono, count = temp2[0], int(temp2[1])
                            mono_count[iupac_comp_mono] = count
                    if len(mono_count) > 0:
                        temp3 = ""
                        for iupac_mono_str, mono_count_eatch in sorted(mono_count.items()):
                            temp3 += iupac_mono_str + str(mono_count_eatch)
                        self.monosaccharide_count[nodeacc] = mono_count

                    for word in l.split()[6:]:
                        m = re.search("^([A-Z]+)(:(.*))?$",word)
                        if m:
                             if m.group(3):
                                 self.extra_properties[nodeacc][m.group(1)] = m.group(3)
                             else:
                                 self.extra_properties[nodeacc][m.group(1)] = True
                             

                elif lineInfo == "edge":
                    fromx,tostr = l.split(':')
                    fromx = fromx.strip()
                    assert(fromx not in content["edges"])
                    content["edges"][fromx] = []
                    for to in tostr.split():
                        if to.endswith("*"):
                            content["edges"][fromx].append((to[:-1],"+aldi"))
                        else:
                            content["edges"][fromx].append((to,"subsumed"))

                elif lineInfo == "tree":
                    gtcacc, g_type, asterik, missingrank = list(re.compile(r"(G\d{5}\w{2}\*?) (BaseComposition|Composition|Topology|Saccharide|-)(\*)? (\d{1,5}|-)").findall(l))[0]
                    if not gtcacc.endswith("*"):
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
            for children in raw_data[m]["edges"].values():
                top = top - set(c[0] for c in children if c[1] == 'subsumed')
            top = [ (n,"subsumed") for n in top ]
            self.alledges[m] = top

        self.allnodestype["00000001"] = "glycan"
        self.alledges["00000001"] = [ (m,"subsumed") for m in allmass ]

        self.allinedges = defaultdict(set)
        for pa,chs in self.alledges.items():
            for ch in chs:
                self.allinedges[ch[0]].add((pa,ch[1]))

        return raw_data

    def root(self):
        return "00000001"

    def nodes(self):
        for n in self.allnodestype.keys():
            yield n

    def parents(self, accession, edgetype=GNOmeAPI.EdgeType.Strict):
        for p in self.allinedges.get(accession, []):
            if edgetype == GNOmeAPI.EdgeType.Any:
                yield p[0]
            elif edgetype == GNOmeAPI.EdgeType.Strict and p[1] == "subsumed":
                yield p[0]
            elif edgetype == GNOmeAPI.EdgeType.Other and p[1] != "subsumed":
                yield p[0]

    def children(self, accession, edgetype=GNOmeAPI.EdgeType.Strict):
        for c in self.alledges.get(accession, []):
            if edgetype == GNOmeAPI.EdgeType.Any:
                yield c[0]
            elif edgetype == GNOmeAPI.EdgeType.Strict and c[1] == "subsumed":
                yield c[0]
            elif edgetype == GNOmeAPI.EdgeType.Other and c[1] != "subsumed":
                yield c[0]

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

    def get_archetype(self, accession):
        retval = self.extra_properties.get(accession,{}).get('HASARCH',None)
        if retval:
            return retval.rstrip('*')
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
        return self.allmissingrank.get(accession)

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

    def generateOWL(self, input_file_path, output_file_path, mass_lut_file_path, acc_file_path, version=None, exact_sym=None, specific_sym=None, replacement=None):

        if specific_sym is None:
            specific_sym = {}
        if exact_sym is None:
            exact_sym = []

        if replacement is None:
            replacement = {}
        assert isinstance(replacement, dict)

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

        r = OWLWriter(mass_LUT_file_path=mass_lut_file_path, historical_accession_file_path=acc_file_path, version=version, replacement=replacement)

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
                node_obj.set_archetype(self.get_archetype(n))


                if n in all_syms:
                    for sym0, sym0_type in all_syms[n]:
                        node_obj.addSynonym(sym0, sym0_type)

            nodes.add(mass)

            for n in nodes:
                for c in self.children(n):
                    r.connect(r.getNode(n), r.getNode(c), "strict")

        for n in self.nodes():
            for c in self.children(n,GNOmeAPI.EdgeType.Other):
                r.connect(r.getNode(n), r.getNode(c), "other")

        f = open(output_file_path, "wb")
        s = r.write(f, r.make_graph())
        f.close()




from rdflib import URIRef, Namespace


class OWLWriter():
    _nodes = {}

    def __init__(self, mass_LUT_file_path=None, historical_accession_file_path=None, version=None, replacement=None):
        self.version = version

        if mass_LUT_file_path:
            self.mass_LUT_file_path = mass_LUT_file_path
        else:
            self.mass_LUT_file_path = "./mass_lut.txt"

        if historical_accession_file_path:
            self.historical_accession_file_path = historical_accession_file_path
        else:
            self.historical_accession_file_path = "./all_accessions.txt"

        if replacement is None:
            self.replacement = {}
        else:
            assert isinstance(replacement, dict)
            self.replacement = replacement

        self.readmassidmap(mass_LUT_file_path)
        self.readaccessionlist(historical_accession_file_path)


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
                x = i.split()
                d[x[1]] = x[0]
        else:
            d["1.01"] = "10000001"
        self.massiddict = d
        self.used_mass = set()

    newMass = False

    def overwritemasslookuptable(self):
        print("new mass ID was assigned")
        mass_lut_file_handle = open(self.mass_LUT_file_path, "w")
        mass_lut_file_handle.write("id\tmass\n")
        for mass in sorted(self.massiddict.keys(),key=float):
            id = self.massiddict[mass]
            mass_lut_file_handle.write("%s\t%s\n" % (id, mass))
        mass_lut_file_handle.close()


    def readaccessionlist(self, fp):
        d = set()

        if fp:
            d = set(open(fp).read().strip().split())

        self.gtc_accession = d
        self.used_gtcacc = set()

    newGTCAcc = False

    def overwriteaccessionlist(self):
        print("new GlyTouCan accession was added")
        file_handle = open(self.historical_accession_file_path, "w")
        for acc in sorted(list(self.gtc_accession)):
            file_handle.write(acc + '\n')
        file_handle.close()


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
    is_subsumed_by_node_annotation_property = 24
    is_restriction_member_annotation_property = 25

    structure_browser_link = 41
    composition_browser_link = 42

    cbbutton_annotation_property = 101
    has_structure_characterization_score = 102

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
        },

        "archetype": {
            "id": 36,
            "label": "has_archetype",
            "definition": """A metadata relation between a glycan and its archetype that is equivalent to the glycan after the removal of reducing end ring information and anomeric configuration."""
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

    def mass_decimal_to_mass_id(self, mass):
        rmass = "%.2f"%(float(mass),)
        try:
            res = self.massiddict[rmass]
        except KeyError:
            res = str(int(max(self.massiddict.values())) + 1)
            self.massiddict[rmass] = res
            self.newMass = True
        return res

    def add_massnode_to_graph(self, graph, node, mass, definition, rdfs, Literal, has_subsumption_level_node, subsumption_level_node, is_subsumed_by_property):
        graph.add((node, rdfs.subClassOf, self.gnouri(self.glycan_class)))
        graph.add((node, is_subsumed_by_property, self.gnouri(self.glycan_class)))
        graph.add((node, has_subsumption_level_node, subsumption_level_node))
        graph.add((node, definition, Literal("A glycan characterized by underivitized molecular weight of %s Daltons." % mass)))
        graph.add((node, rdfs.label, Literal("glycan of molecular weight %s Da" % mass)))

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

        is_subsumed_by_node = self.gnouri(self.is_subsumed_by_node_annotation_property)

        outputGraph.add((is_subsumed_by_node, rdf.type, owl.AnnotationProperty))
        outputGraph.add((is_subsumed_by_node, rdfs.label, Literal("is_subsumed_by")))
        outputGraph.add((is_subsumed_by_node, definition,
                         Literal("A metadata relation between a glycan and a glycan that subsumes it. Note that this relation is a superset of the subClassOf relation as molecular weight is not necessarily conserved as it is for subClassOf.")))

        is_restriction_member = self.gnouri(self.is_restriction_member_annotation_property)

        outputGraph.add((is_restriction_member, rdf.type, owl.AnnotationProperty))
        outputGraph.add((is_restriction_member, rdfs.label, Literal("is_restriction_member")))
        outputGraph.add((is_restriction_member, definition,
                         Literal("Whether the glycan is a member of the restriction. If not present, the glycan is not considered to be in the restriction.")))

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
        for sl in ("basecomposition", "composition", "topology", "archetype"):
            has_xxx_node = self.gnouri(self.has_subsumption_level[sl]["id"])

            outputGraph.add((has_xxx_node, rdf.type, owl.AnnotationProperty))
            outputGraph.add((has_xxx_node, rdfs.label, Literal(self.has_subsumption_level[sl]["label"])))
            outputGraph.add((has_xxx_node, definition, Literal(self.has_subsumption_level[sl]["definition"])))

            has_xxx_nodes[sl] = has_xxx_node

        has_topology_node = has_xxx_nodes["topology"]
        has_composition_node = has_xxx_nodes["composition"]
        has_basecomposition_node = has_xxx_nodes["basecomposition"]
        has_archetype_node = has_xxx_nodes["archetype"]

        # Add AnnotationProperty for IUPAC composition (for the viewer)
        cbbutton_node = self.gnouri(self.cbbutton_annotation_property)

        outputGraph.add((cbbutton_node, rdf.type, owl.AnnotationProperty))
        outputGraph.add((cbbutton_node, rdfs.label, Literal("_widget_button_state")))
        outputGraph.add((cbbutton_node, definition, Literal("Internal raw data for both GNOme browser.")))

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

        # Add AnnotationProperty for has_structure_characterization_score
        has_structure_characterization_score_node = self.gnouri(self.has_structure_characterization_score)

        outputGraph.add((has_structure_characterization_score_node, rdf.type, owl.AnnotationProperty))
        outputGraph.add((has_structure_characterization_score_node, definition, Literal(
            "A score for the extent of characterization provided by the glycan's description. Glycan descriptions that completely characterize a glycan have score 0. Scores increase monotonically with subsumption. Scores should only be compared between glycan descriptions with the same monosaccharide (base-)composition. Score may change in future releases.")))
        outputGraph.add((has_structure_characterization_score_node, rdfs.label, Literal("has_structure_characterization_score")))


        # Add AnnotationProperty for consider
        consider_node = oboInOwl["consider"]

        outputGraph.add((consider_node, rdf.type, owl.AnnotationProperty))
        outputGraph.add((consider_node, rdfs.isDefinedBy, oboInOwl[""]))
        outputGraph.add((consider_node, rdfs.label, Literal("consider")))

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
                rdfNode = self.gnouri(self.mass_decimal_to_mass_id(n.getID()))

            outputGraph.add((rdfNode, rdf.type, owl.Class))

            if not n.ismolecularweight():
                self.used_gtcacc.add(n.getID())

                outputGraph.add((rdfNode, has_glytoucan_id_node, Literal(n.getID())))
                outputGraph.add((rdfNode, has_glytoucan_link_node, gtcs[n.getID()]))
                outputGraph.add((rdfNode, definition,
                                 Literal("A glycan described by the GlyTouCan entry with accession %s." % n.getID())))
                outputGraph.add((rdfNode, rdfs.label,
                                 Literal("%s" % n.getID())))

                has_structure_characterization_score_node
                outputGraph.add((rdfNode, has_structure_characterization_score_node, Literal(n.missing_rank())))

                assert n._nodeType
                outputGraph.add((rdfNode, has_subsumption_level_node, subsumptionLevel[n._nodeType.lower()]))

            else:
                self.used_mass.add(n.getID())
                self.add_massnode_to_graph(outputGraph, rdfNode, n.getID(), definition, rdfs, Literal, has_subsumption_level_node, subsumptionLevel["molecularweight"], is_subsumed_by_node)

            for sym, sym_type in n.synonym():
                outputGraph.add((rdfNode, sym_types[sym_type], Literal(sym)))

            bcomp_acc = n.get_basecomposition()
            comp_acc = n.get_composition()
            topo_acc = n.get_topology()
            arch_acc = n.get_archetype()
            if bcomp_acc:
                outputGraph.add((rdfNode, has_basecomposition_node, self.gnouri(bcomp_acc)))
            if comp_acc:
                outputGraph.add((rdfNode, has_composition_node, self.gnouri(comp_acc)))
            if topo_acc:
                outputGraph.add((rdfNode, has_topology_node, self.gnouri(topo_acc)))
            if arch_acc:
                outputGraph.add((rdfNode, has_archetype_node, self.gnouri(arch_acc)))

            cbbutton_str = n.get_iupac_composition()
            if cbbutton_str:
                outputGraph.add((rdfNode, cbbutton_node, Literal(cbbutton_str)))

            if not n.ismolecularweight():
                if n.isbasecomposition() or n.iscomposition():
                    outputGraph.add((rdfNode, has_composition_browser_node, URIRef(
                        "https://gnome.glyomics.org/CompositionBrowser.html?focus=%s" % n.getID())))
                outputGraph.add((rdfNode, has_structure_browser_node, URIRef(
                    "https://gnome.glyomics.org/StructureBrowser.html?focus=%s" % n.getID())))


        for l in self.allRelationship():
            if l.getEdgeType() in ("strict","other"):
                if l._sideA._nodeType == "molecularweight":
                    id = self.massiddict["%.2f"%(float(l._sideA.getID()),)]
                    n1 = self.gnouri(id)
                else:
                    n1 = self.gnouri(l._sideA.getID())
                n2 = self.gnouri(l._sideB.getID())
                if l.getEdgeType() == "strict":
                    outputGraph.add((n2, rdfs.subClassOf, n1))
                    outputGraph.add((n2, is_subsumed_by_node, n1))
                elif l.getEdgeType() == "other":
                    outputGraph.add((n2, is_subsumed_by_node, n1))

        tmp1 = len(self.gtc_accession)
        self.gtc_accession = self.gtc_accession.union(self.used_gtcacc)
        tmp2 = len(self.gtc_accession)
        if tmp2 > tmp1:
            self.newGTCAcc = True

        if self.newMass:
            self.overwritemasslookuptable()

        if self.newGTCAcc:
            self.overwriteaccessionlist()



        unused_mass = set(self.massiddict.keys()) - set(map(lambda m: "%.2f"%(float(m),), self.used_mass))
        unused_gtc_acc = self.gtc_accession - self.used_gtcacc

        rdfNodeXSDTrue = rdflib.Literal("true", datatype=rdflib.XSD.boolean)
        for n in unused_mass:
            rdfNode = self.gnouri(self.mass_decimal_to_mass_id(n))
            outputGraph.add((rdfNode, rdf.type, owl.Class))
            self.add_massnode_to_graph(outputGraph, rdfNode, n, definition, rdfs, Literal, has_subsumption_level_node, subsumptionLevel["molecularweight"],is_subsumed_by_node)

        for n in unused_gtc_acc:
            rdfNode = self.gnouri(n)
            outputGraph.add((rdfNode, rdf.type, owl.Class))
            outputGraph.add((rdfNode, owl.deprecated, rdfNodeXSDTrue))

            outputGraph.add((rdfNode, definition,
                             Literal("OBSOLETE. A glycan described by the GlyTouCan entry with accession %s." % n)))
            outputGraph.add((rdfNode, rdfs.label,
                             Literal("obsolete %s" % n)))



        for oldacc, newacc in self.replacement.items():
            if oldacc in self.gtc_accession and newacc in self.used_gtcacc:
                rdfNodeStringNewAcc = rdflib.Literal(newacc, datatype=rdflib.XSD.string)
                outputGraph.add((self.gnouri(oldacc), consider_node, rdfNodeStringNewAcc))

        return outputGraph

    def write(self, handle, graph):
        from rdflib.plugins.serializers.rdfxml import PrettyXMLSerializer
        writer = PrettyXMLSerializer(SubjectOrderedGraph(graph), max_depth=1)
        writer.store.base = None
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

    def set_archetype(self, archetype):
        self._archetype = archetype

    def get_archetype(self):
        try:
            return self._archetype
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
    _edgeTypes = tuple(["equal", "strict", "other"])
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

    def get_accessions(self, *restriction_set_names):
        accs = dict()
        for n in restriction_set_names:
            url = self.restriction_url % n
            accs.update(dict(sorted(map(lambda l: (l[0],l[-1]),map(lambda l: l.split(),open(url).readlines())))))
        n1 = len(accs)
        n2 = sum(1 for acc in accs if accs[acc] == acc)
        if n1 == n2:
            return list(sorted(accs.keys()))
        return dict(sorted(accs.items()))

    def getoutputpath(self):
        raise NotImplemented()

    def write(self, output_path):
        d = self.getdata()
        json.dump(d, open(output_path, "w"), indent=2, separators=(',',': '), sort_keys=True)


class GNOme_Theme_GlyGen(GNOme_Theme_Base):

    def getdata(self):
        return {
            "icon_style": "snfg",
            "image_url_prefix": "https://glymage.glyomics.org/image/snfg/extended/",
            "image_url_suffix": ".svg",
            "brand": "A <a href='https://www.glygen.org/' target='_blank'>GlyGen</a> project",
            "external_resources": [
                {
                    "name": "GlyGen",
                    "url_prefix": "https://www.glygen.org/glycan/",
                    "url_suffix": "",
                    "glycan_set": self.get_accessions("GlyGen")
                }
            ]

        }

class GNOme_Theme_PubChemCID(GNOme_Theme_Base):

    def getdata(self):
        return {
            "icon_style": "snfg",
            "image_url_prefix": "https://glymage.glyomics.org/image/snfg/extended/",
            "image_url_suffix": ".svg",
            "external_resources": [
                {
                    "name": "PubChem",
                    "url_prefix": "https://pubchem.ncbi.nlm.nih.gov/compound/",
                    "url_suffix": "",
                    "glycan_set": self.get_accessions("PubChemCID")
                },{
                    "name": "GlyTouCan",
                    "url_prefix": "https://glytoucan.org/Structures/Glycans/",
                    "url_suffix": "",
                    "glycan_set": None
                }
            ]

        }

class GNOme_Theme_GlyConnect(GNOme_Theme_Base):

    def getdata(self):
        return {
            "icon_style": "snfg",
            "image_url_prefix": "https://glymage.glyomics.org/image/snfg/extended/",
            "image_url_suffix": ".svg",
            "external_resources": [
                {
                    "name": "GlyConnect",
                    "url_prefix": "https://glyconnect.expasy.org/all/structures/",
                    "url_suffix": "",
                    "glycan_set": self.get_accessions("GlyConnect")
                },{
                    "name": "GlyTouCan",
                    "url_prefix": "https://glytoucan.org/Structures/Glycans/",
                    "url_suffix": "",
                    "glycan_set": None
                }
            ]

        }

class GNOme_Theme_GlyGen_Sandbox(GNOme_Theme_Base):

    def getdata(self):
        return {
            "icon_style": "snfg",
            "image_url_prefix": "https://glymage.glyomics.org/image/snfg/extended/",
            "image_url_suffix": ".svg",
            "brand": "A <a href='https://www.glygen.org/' target='_blank'>GlyGen</a> project",
            "external_resources": [
                {
                    "name": "GlyGen",
                    "url_prefix": "https://www.glygen.org/glycan/",
                    "url_suffix": "",
                    "glycan_set": self.get_accessions("GlyGen")
                },
                {
                    "name": "Sandbox",
                    "url_prefix": "https://sandbox.glyomics.org/explore.html?focus=",
                    "url_suffix": "",
                    "glycan_set": self.get_accessions("GlycoTree_NGlycans","GlycoTree_OGlycans")
                }
            ]

        }

class GNOme_Theme_Default(GNOme_Theme_Base):

    def getdata(self):
        return {
            "icon_style": "snfg",
            "image_url_prefix": "https://glymage.glyomics.org/image/snfg/extended/",
            "image_url_suffix": ".svg",
            "brand": None,
            "external_resources": [
                {
                    "name": "GlyGen",
                    "url_prefix": "https://www.glygen.org/glycan/",
                    "url_suffix": "",
                    "glycan_set": self.get_accessions("GlyGen")
                },{
                    "name": "GlyCosmos",
                    "url_prefix": "https://glycosmos.org/glycans/show/",
                    "url_suffix": "",
                    "glycan_set": self.get_accessions("GlyCosmos")
                },{
                    "name": "GlyTouCan",
                    "url_prefix": "https://glytoucan.org/Structures/Glycans/",
                    "url_suffix": "",
                    "glycan_set": None
                }
            ]

        }


class IncompleteScore:

    def __init__(self):
        pass

    def monoscore(self, m, subst_list):
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
            res += (float(config.count(None)) / len(config)) * 2
        else:
            res += 2

        if anomer != Anomer.uncyclized:
            if rs == None:
                res += 1
            if re == None:
                res += 1

        if stem == None or None in stem:
            res += 2

        #if superclass == None:
        #    res += 2

        subst_score_total = 0.
        for sl in m.substituent_links():
            subst = sl.child()

            if subst.name() == Substituent.nAcetyl or\
                    (subst.name() == Substituent.nglycolyl and sl.parent().superclass() == 9):

                subst_score_total += 2
                if sl.parent_pos() == None:
                    res += 2
                else:
                    l = len(sl.parent_pos())
                    if l > 1:
                        res += 1

                subst_list.remove(subst)


        res = res / 8 * 10
        return res

    def linksimplescore(self, l):
        res = 0.0
        pp = l.parent_pos()
        cp = l.child_pos()
        pt = l.parent_type()
        ct = l.child_type()

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

        res = res / 8. * 10 # scale from 8 to 10
        return res


    def substscore(self, s):
        if s._sub != None:
            return 2
        return 0

    def score(self, g):

        monos = list(g.all_nodes())
        total_mono = len(monos)
        unconnected_mono_root = list(g.unconnected_roots())
        unconnected_mono_root = list(filter(lambda m: m.is_monosaccharide(), unconnected_mono_root))
        det_parent_links = [list(m.parent_links()) for m in filter(lambda m: m not in unconnected_mono_root and m != g.root(), monos)]
        und_parent_links = [list(m.parent_links()) for m in unconnected_mono_root]

        allsubst = list(filter(lambda n: not n.is_monosaccharide(), g.all_nodes(subst=True, undet_subst=False)))

        mono_score_result = 0.0
        mono_score_total  = 10.0 * total_mono

        link_score_result = 0.0
        link_score_total  = 10.0 * (total_mono - 1)

        for m in monos:
            mono_score_result += self.monoscore(m, allsubst)

        linkage_simple_parameter_weight = 0.65
        linkage_undetermined_weight = 1 - linkage_simple_parameter_weight
        if len(und_parent_links) == total_mono:
            # Basecomp...
            link_score_result = 10.0 * (total_mono - 1)
        else:
            for ls in und_parent_links + det_parent_links:

                if len(ls) == 1:
                    link_score_result += self.linksimplescore(ls[0]) * linkage_simple_parameter_weight

                elif len(ls) > 1:
                    parent_link_scores = []
                    for pl in ls:
                        s = self.linksimplescore(pl) * linkage_simple_parameter_weight
                        parent_link_scores.append(s)

                    ave = sum(parent_link_scores) / len(parent_link_scores)
                    und_ratio = (float(len(parent_link_scores)-1) / (total_mono))
                    s = ave + 10*linkage_undetermined_weight * und_ratio
                    link_score_result += s

                else:
                    # Unknown of from-and-to and detail of the link
                    link_score_result += 10.0




        for subst in allsubst:
            if str(subst) == "anhydro":
                continue

            pl = list(subst.parent_links())
            link_score_total += 5

            if len(pl) == 1:
                link_score_result += self.linksimplescore(pl[0]) * linkage_simple_parameter_weight/2

            elif len(pl) > 1:

                parent_link_scores = []
                for pl0 in pl:
                    s = self.linksimplescore(pl0) * linkage_simple_parameter_weight/2
                    parent_link_scores.append(s)

                ave = sum(parent_link_scores) / len(parent_link_scores)
                und_ratio = (float(len(parent_link_scores)-1) / (total_mono))
                s = ave + 5*linkage_undetermined_weight * und_ratio
                link_score_result += s

            else:
                link_score_result += 5.0


        if link_score_total == 0:
            if link_score_result !=  0:
                raise RuntimeError
            link_score_total = 1


        monosaccharide_weight = 1 / (1 + 4. / (5 * linkage_simple_parameter_weight))
        monosaccharide_weight = linkage_simple_parameter_weight / (1 + linkage_simple_parameter_weight)
        linkage_weight = 1 - monosaccharide_weight

        #
        #v1 = 2./8 * monosaccharide_weight
        #v2 = 2./8 * linkage_simple_parameter_weight * linkage_weight
        #print monosaccharide_weight, linkage_weight, v1, v2

        res = mono_score_result/mono_score_total * monosaccharide_weight + link_score_result/link_score_total * linkage_weight
        #print mono_score_result / mono_score_total, link_score_result / link_score_total
        res = int(res * 10000)
        return res

def main():

    cmd = sys.argv[1]
    sys.argv.pop(1)

    if cmd == "restrict":

        restriction = set(map(lambda l: l.split(None,1)[0],open(sys.argv[1]).read().splitlines()))

        g = GNOme()
        g.restrict(restriction)
        g.write(sys.stdout)

    elif cmd == "dump":

        g = GNOme()
        g.dump()

    elif cmd == "compute":

        verbose = 0
        usecache=False
        while len(sys.argv) > 1 and sys.argv[1] in ("-v","-c"):
            if sys.argv[1] == "-v":
                verbose += 1
                sys.argv.pop(1)
            elif sys.argv[1] == "-c":
                usecache=True
                sys.argv.pop(1)

        g = SubsumptionGraph()
        g.compute(*sys.argv[1:], verbose=verbose, usecache=usecache)

    elif cmd == "writeowl":

        kv_para = {
            "version": None,
            "archive": None,
            "replace": None
        }
        if len(sys.argv) < 5:
            print("Please provide dumpfile, output file path(with file name), mass LUT path and version (optional)")
            sys.exit(1)
        if len(sys.argv) > 5:
            for k, v in zip(sys.argv[5::2], sys.argv[6::2]):
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

        replace_mapping = {}

        if kv_para["archive"] is not None:
            archive_file_handle = open(kv_para["archive"])
            for i, l in enumerate(archive_file_handle):
                acc = l.strip()
                if i == 0 and "accession" in acc:
                    continue
                replace_mapping[acc] = None


        if kv_para["replace"] is not None:
            replace_file_path = kv_para["replace"]

            replace_file_handle = open(replace_file_path)
            for i, l in enumerate(replace_file_handle):
                linfo = l.strip().split()
                if i == 0 and "accession" in linfo:
                    continue

                newacc, retire = linfo
                replace_mapping[retire] = newacc



        #for syms in exactSym:
        #    print len(syms)
        #for k, syms in specificSym.items():
        #    print k, len(syms)

        ifn = sys.argv[1]  # path to input file
        ofn = sys.argv[2]
        mass_lut = sys.argv[3]
        accession_list_file = sys.argv[4]

        subsumption_instance = SubsumptionGraph()
        subsumption_instance.generateOWL(ifn, ofn, mass_lut, accession_list_file, version=kv_para["version"], exact_sym=exactSym, specific_sym=specificSym, replacement=replace_mapping)

        if "allExactSymOutput" in kv_para:
            allexactsym = {}
            for symset in exactSym + list(specificSym.values()):
                for acc, sym0 in symset.items():
                    if acc not in allexactsym:
                        allexactsym[acc] = []
                    allexactsym[acc].append(sym0)

            allExactSymOutputF = open(kv_para["allExactSymOutput"], "w")
            for acc in sorted(allexactsym.keys()):
                for sym0 in allexactsym[acc]:
                    allExactSymOutputF.write("%s\t%s\n" % (acc, sym0))
            allExactSymOutputF.close()


    elif cmd in ("writeresowl","writeresowl_with_ancestor_structures"):

        if len(sys.argv) < 4:
            print("Please provide GNOme.owl, restriction set name, output file path")
            sys.exit(1)

        ifn = sys.argv[1]
        ofn = sys.argv[3]
        restriction_accs_file = sys.argv[2]

        accs = set(map(lambda l: l.split(None,1)[0],open(restriction_accs_file).read().splitlines()))

        GNOme_res = GNOme(resource=ifn)
        if cmd == "writeresowl":
            GNOme_res.restrict(accs)
        elif cmd == "writeresowl_with_ancestor_structures":
            GNOme_res.restrict(accs,allstructanc=True,alltopoanc=True)

        f = open(ofn, "wb")
        GNOme_res.write(f)
        f.close()

    elif cmd == "viewerdata":
        # python GNOme.py viewerdata ./GNOme.owl ./gnome_subsumption_raw.txt ./GNOme.browser.js

        if len(sys.argv) < 3:
            print("Please provide GNOme.owl and output file path")
            sys.exit(1)

        ifn = sys.argv[1]
        ofn_data = sys.argv[2]
        gnome = GNOme(resource=ifn)
        gnome.toViewerData(ofn_data)


    elif cmd == "UpdateAcc":

        if len(sys.argv) < 4:
            print("Please provide restriction set name and file path")
            sys.exit(1)

        restriction_set_name = sys.argv[1]
        fp = sys.argv[2]

        restriction_set = []
        if restriction_set_name == "BCSDB":
            sys.exit(0)

        #elif restriction_set_name == "GlyGen":
        #    fp = os.path.dirname(os.path.abspath(__file__)) + "/../smw/glycandata/data/glygen_accessions.txt"
        #    restriction_set = open(fp).read().strip().split()

        elif restriction_set_name in ["GlyGen"]:
            glycandata_tsv_fp = "../smw/glycandata/export/allglycan.tsv"

            restriction_set = open(glycandata_tsv_fp).read().strip().split()
            restriction_set.pop(0)
            restriction_set = sorted(restriction_set)
        else:
            print("Restriction set: %s is not supported")
            sys.exit(1)

        open(fp, "w").write("\n".join(restriction_set))

        json_fp = open(sys.argv[3], "w")
        json.dump(restriction_set, json_fp, indent=2, separators=(',',': '), sort_keys=True)

    elif cmd == "UpdateTheme":

        if len(sys.argv) < 3:
            print("Please provide restriction set name and file path")
            sys.exit(1)

        restriction_url = sys.argv[1]
        theme_path = sys.argv[2]

        td = GNOme_Theme_Default(restriction_url)
        td.write(theme_path + "default.json")

        tgg = GNOme_Theme_GlyGen(restriction_url)
        tgg.write(theme_path + "GlyGen.json")

        tggs = GNOme_Theme_GlyGen_Sandbox(restriction_url)
        tggs.write(theme_path + "Sandbox.json")

        tpc = GNOme_Theme_PubChemCID(restriction_url)
        tpc.write(theme_path + "PubChemCID.json")

        tpc = GNOme_Theme_GlyConnect(restriction_url)
        tpc.write(theme_path + "GlyConnect.json")

    else:

        gnome = GNOme()

        if hasattr(gnome,cmd):
            result = getattr(gnome,cmd)(*sys.argv[1:])
            if isinstance(result,basestring):
                print(result)
            else:
                for r in result:
                    print(r)
        else:
            print >> sys.stderr, "Bad command: %s" % (cmd,)
            sys.exit(1)

if __name__ == "__main__":
    main()
