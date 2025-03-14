#!/bin/env python3.12

import sys
from collections import defaultdict

from getwiki import GlycanData, Glycan
w = GlycanData()

# print >>sys.stderr, "Read subsumption graph"
from pygly.GNOme import SubsumptionGraph

gnome = SubsumptionGraph()
gnome.loaddata(sys.argv[1])

gnomeacc = set()
for n in gnome.nodes():
    gnomeacc.add(n)

accs = set()
for fn in sys.argv[2:]:
    accs.update(open(fn).read().split())

for m in w.iterglycan():

    acc = m.get('accession')
    m.delete_annotations(type="Subsumption")

    if acc in gnomeacc:
        m.set_annotation(value=acc,property="GNOme",source="GNOme",type="CrossReference")
    else:
        m.delete_annotations(property="GNOme",source="GNOme",type="CrossReference")
 
    if gnome.isbasecomposition(acc):
        m.set_annotation(property="Level",value="BaseComposition",source="GNOme",type="Subsumption")
        m.set_annotation(property="BaseComposition",value=acc,source="GNOme",type="Subsumption")
        hasbcomp = list(gnome.has_basecomposition(acc))
        comps = list(filter(accs.__contains__,filter(gnome.iscomposition,hasbcomp)))
        m.set_annotation(property="Compositions",value=comps,source="GNOme",type="Subsumption")
        topos = list(filter(accs.__contains__,filter(gnome.istopology,hasbcomp)))
        m.set_annotation(property="Topologies",value=topos,source="GNOme",type="Subsumption")
        saccs = list(filter(accs.__contains__,filter(gnome.issaccharide,hasbcomp)))
        m.set_annotation(property="Saccharides",value=saccs,source="GNOme",type="Subsumption")
    elif gnome.iscomposition(acc):
        m.set_annotation(property="Level",value="Composition",source="GNOme",type="Subsumption")
        bcomp = gnome.get_basecomposition(acc)
        m.set_annotation(property="BaseComposition",value=bcomp,source="GNOme",type="Subsumption")
        m.set_annotation(property="Composition",value=acc,source="GNOme",type="Subsumption")
        hascomp = list(gnome.has_composition(acc))
        topos = list(filter(gnome.istopology,hascomp))
        m.set_annotation(property="Topologies",value=topos,source="GNOme",type="Subsumption")
        saccs = list(filter(gnome.issaccharide,hascomp))
        m.set_annotation(property="Saccharides",value=saccs,source="GNOme",type="Subsumption")
    elif gnome.istopology(acc):
        m.set_annotation(property="Level",value="Topology",source="GNOme",type="Subsumption")
        bcomp = gnome.get_basecomposition(acc)
        m.set_annotation(property="BaseComposition",value=bcomp,source="GNOme",type="Subsumption")
        comp = gnome.get_composition(acc)
        m.set_annotation(property="Composition",value=comp,source="GNOme",type="Subsumption")
        m.set_annotation(property="Topology",value=acc,source="GNOme",type="Subsumption")
        hastopo = list(gnome.has_topology(acc))
        saccs = list(filter(gnome.issaccharide,hastopo))
        m.set_annotation(property="Saccharides",value=saccs,source="GNOme",type="Subsumption")
    elif gnome.issaccharide(acc):
        m.set_annotation(property="Level",value="Saccharide",source="GNOme",type="Subsumption")
        bcomp = gnome.get_basecomposition(acc)
        m.set_annotation(property="BaseComposition",value=bcomp,source="GNOme",type="Subsumption")
        comp = gnome.get_composition(acc)
        m.set_annotation(property="Composition",value=comp,source="GNOme",type="Subsumption")
        topo = gnome.get_topology(acc)
        m.set_annotation(property="Topology",value=topo,source="GNOme",type="Subsumption")

    m.set_annotation(property="Subsumes",value=list(filter(accs.__contains__,gnome.children(acc))),source="GNOme",type="Subsumption")
    m.set_annotation(property="SubsumedBy",value=list(filter(accs.__contains__,gnome.parents(acc))),source="GNOme",type="Subsumption")
    m.set_annotation(property="Descendant",value=list(filter(accs.__contains__,gnome.descendants(acc))),source="GNOme",type="Subsumption")
    m.set_annotation(property="Ancestor",value=list(filter(accs.__contains__,gnome.ancestors(acc))),source="GNOme",type="Subsumption")
    m.set_annotation(property="MissingScore",value=gnome.get_missing_rank(acc),source="GNOme",type="Subsumption")
    m.set_annotation(property="FullyDetermined",value=list(filter(lambda acc1: acc1 in accs and int(gnome.get_missing_rank(acc1)) == 0,gnome.descendants(acc))),source="GNOme",type="Subsumption")
    m.set_annotation(property="Leaf",value=list(filter(lambda acc1: acc1 in accs and gnome.isleaf(acc1),gnome.descendants(acc))),source="GNOme",type="Subsumption")
    
    if w.put(m):
        print(acc)
