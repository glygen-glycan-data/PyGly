#!/bin/env python27

import sys
from collections import defaultdict

from getwiki import GlycanData, Glycan
w = GlycanData()

# print >>sys.stderr, "Read subsumption graph"
from pygly.GNOme import SubsumptionGraph

gnome = SubsumptionGraph()
gnome.loaddata(sys.argv[1])

accs = set()
for fn in sys.argv[2:]:
    accs.update(open(fn).read().split())

for m in w.iterglycan():

    acc = m.get('accession')
    m.delete_annotations(type="Subsumption")

    if gnome.isbasecomposition(acc):
	m.set_annotation(property="Level",value="BaseComposition",source="GNOme",type="Subsumption")
	m.set_annotation(property="BaseComposition",value=acc,source="GNOme",type="Subsumption")
	hasbcomp = list(gnome.has_basecomposition(acc))
        comps = filter(accs.__contains__,filter(gnome.iscomposition,hasbcomp))
	m.set_annotation(property="Compositions",value=comps,source="GNOme",type="Subsumption")
        topos = filter(accs.__contains__,filter(gnome.istopology,hasbcomp))
	m.set_annotation(property="Topologies",value=topos,source="GNOme",type="Subsumption")
        saccs = filter(accs.__contains__,filter(gnome.issaccharide,hasbcomp))
	m.set_annotation(property="Saccharides",value=saccs,source="GNOme",type="Subsumption")
    elif gnome.iscomposition(acc):
	m.set_annotation(property="Level",value="Composition",source="GNOme",type="Subsumption")
	bcomp = gnome.get_basecomposition(acc)
	m.set_annotation(property="BaseComposition",value=bcomp,source="GNOme",type="Subsumption")
	m.set_annotation(property="Composition",value=acc,source="GNOme",type="Subsumption")
	hascomp = list(gnome.has_composition(acc))
        topos = filter(gnome.istopology,hascomp)
	m.set_annotation(property="Topologies",value=topos,source="GNOme",type="Subsumption")
        saccs = filter(gnome.issaccharide,hascomp)
	m.set_annotation(property="Saccharides",value=saccs,source="GNOme",type="Subsumption")
    elif gnome.istopology(acc):
	m.set_annotation(property="Level",value="Topology",source="GNOme",type="Subsumption")
	bcomp = gnome.get_basecomposition(acc)
	m.set_annotation(property="BaseComposition",value=bcomp,source="GNOme",type="Subsumption")
	comp = gnome.get_composition(acc)
	m.set_annotation(property="Composition",value=comp,source="GNOme",type="Subsumption")
	m.set_annotation(property="Topology",value=acc,source="GNOme",type="Subsumption")
	hastopo = list(gnome.has_topology(acc))
        saccs = filter(gnome.issaccharide,hastopo)
	m.set_annotation(property="Saccharides",value=saccs,source="GNOme",type="Subsumption")
    elif gnome.issaccharide(acc):
	m.set_annotation(property="Level",value="Saccharide",source="GNOme",type="Subsumption")
	bcomp = gnome.get_basecomposition(acc)
	m.set_annotation(property="BaseComposition",value=bcomp,source="GNOme",type="Subsumption")
	comp = gnome.get_composition(acc)
	m.set_annotation(property="Composition",value=comp,source="GNOme",type="Subsumption")
	topo = gnome.get_topology(acc)
	m.set_annotation(property="Topology",value=topo,source="GNOme",type="Subsumption")

    m.set_annotation(property="Subsumes",value=filter(accs.__contains__,gnome.children(acc)),source="GNOme",type="Subsumption")
    m.set_annotation(property="SubsumedBy",value=filter(accs.__contains__,gnome.parents(acc)),source="GNOme",type="Subsumption")
    m.set_annotation(property="Descendant",value=filter(accs.__contains__,gnome.descendants(acc)),source="GNOme",type="Subsumption")
    m.set_annotation(property="Ancestor",value=filter(accs.__contains__,gnome.ancestors(acc)),source="GNOme",type="Subsumption")
    
    if w.put(m):
        print acc
