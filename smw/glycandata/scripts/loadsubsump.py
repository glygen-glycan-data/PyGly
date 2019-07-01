#!/bin/env python27

import sys
from collections import defaultdict

from getwiki import GlycanData, Glycan
w = GlycanData()

# print >>sys.stderr, "Read subsumption graph"

cat = dict()
comp = dict()
topo = dict()
bcomp = dict()
subsumes = defaultdict(set)
sec = None
for l in open(sys.argv[1]):
    if l.startswith('#'):
	sec = l.split()[1]
	continue
    if sec == "NODES":
	sl = l.split()
	acc = sl[0]
	if sl[3] != "None":
            topo[acc] = (sl[3].strip("*"),sl[3].endswith("*"))
	else:
	    topo[acc] = (None,False)
	if sl[4] != "None":
            comp[acc] = (sl[4].strip("*"),sl[4].endswith("*"))
        else:
	    comp[acc] = (None,False)
	if sl[5] != "None":
            bcomp[acc] = (sl[5].strip("*"),sl[5].endswith("*"))
	else:
	    bcomp[acc] = (None,False)
        cat[acc] = (sl[2].strip("*"),sl[2].endswith("*"))
	continue
    if sec == "EDGES":
        sl = l.split()
        fr = sl[0][:-1]
        tos = sl[1:]
        subsumes[fr].update(tos)
        
# print >>sys.stderr, "Inverse subsumption groups"

hascomp = defaultdict(set)
for acc,(c,edlab) in comp.iteritems():
    hascomp[c].add(acc)
hastopo = defaultdict(set)
for acc,(t,edlab) in topo.iteritems():
    hastopo[t].add(acc)
hasbcomp = defaultdict(set)
for acc,(bc,edlab) in bcomp.iteritems():
    hasbcomp[bc].add(acc)

# print >>sys.stderr, "Read restriction(s)"

restriction = None
for f in sys.argv[2:]:
    if not restriction:
        restriction = set()
    restriction.update(open(f).read().split())

def alldescendents(acc):
    desc = set()
    for to in subsumes[acc]:
        if to not in desc:
            desc.add(to)
            desc.update(alldescendents(to))
    return desc

if restriction != None:

    # print >>sys.stderr, "Add all descendents"

    for fr in list(subsumes):
        subsumes[fr].update(alldescendents(fr))

    # print >>sys.stderr, "Restrict subsumptions"

    for fr in list(subsumes):
        if fr not in restriction:
            del subsumes[fr]
        else:
            subsumes[fr] = subsumes[fr]&restriction

    # print >>sys.stderr, "Remove shortcuts"

    toremove = set()
    for n1 in subsumes:
        for n2 in subsumes[n1]:
            for n3 in subsumes[n2]:
                if n3 in subsumes[n1]:
                    toremove.add((n1,n3))

    for n1,n3 in toremove:
        subsumes[n1].remove(n3)

# print >>sys.stderr, "Determine subsumedby relationship(s)"

subsumedby = defaultdict(set)
for fr,tos in subsumes.items():
    for to in tos:
        subsumedby[to].add(fr)

from functools import partial

def keepcat(c,acc):
    return (cat.get(acc)==c and (restriction == None or acc in restriction))

# print >>sys.stderr, "Push to glycandata"

for m in w.iterglycan():

    acc = m.get('accession')
    if cat.get(acc,(None,False))[0] == "BaseComposition":
        comps = filter(partial(keepcat,"Composition"),hasbcomp[acc])
        m.set_annotation(property="Compositions",value=comps,source="EdwardsLab",type="Subsumption")
        topos = filter(partial(keepcat,"Topology"),hasbcomp[acc])
        m.set_annotation(property="Topologies",value=topos,source="EdwardsLab",type="Subsumption")
        saccs = filter(partial(keepcat,"Saccharide"),hasbcomp[acc])
        m.set_annotation(property="Saccharides",value=saccs,source="EdwardsLab",type="Subsumption")
    elif cat.get(acc,(None,False))[0] == "Composition":
        topos = filter(partial(keepcat,"Topology"),hascomp[acc])
        m.set_annotation(property="Topologies",value=topos,source="EdwardsLab",type="Subsumption")
        saccs = filter(partial(keepcat,"Saccharide"),hascomp[acc])
        m.set_annotation(property="Saccharides",value=saccs,source="EdwardsLab",type="Subsumption")
    elif cat.get(acc,(None,False))[0] == "Topology":
        saccs = filter(partial(keepcat,"Saccharide"),hastopo[acc])
        m.set_annotation(property="Saccharides",value=saccs,source="EdwardsLab",type="Subsumption")

    if topo.get(acc,(None,False))[1] and not m.has_annotations(property="Topology",source="GlyTouCan",type="Subsumption"):
	m.set_annotation(property="Topology",value=topo.get(acc)[0],source="EdwardsLab",type="Subsumption")
    else:
	m.delete_annotations(property="Topology",source="EdwardsLab",type="Subsumption")

    if comp.get(acc,(None,False))[1] and not m.has_annotations(property="Composition",source="GlyTouCan",type="Subsumption"):
	m.set_annotation(property="Composition",value=comp.get(acc)[0],source="EdwardsLab",type="Subsumption")
    else:
	m.delete_annotations(property="Composition",source="EdwardsLab",type="Subsumption")

    if bcomp.get(acc,(None,False))[1] and not m.has_annotations(property="BaseComposition",source="GlyTouCan",type="Subsumption"):
	m.set_annotation(property="BaseComposition",value=bcomp.get(acc)[0],source="EdwardsLab",type="Subsumption")
    else:
	m.delete_annotations(property="BaseComposition",source="EdwardsLab",type="Subsumption")

    if cat.get(acc,(None,False))[1] and not m.has_annotations(property="SubsumptionLevel",source="GlyTouCan",type="Subsumption"):
	m.set_annotation(property="SubsumptionLevel",value=cat.get(acc)[0],source="EdwardsLab",type="Subsumption")
    else:
	m.delete_annotations(property="SubsumptionLevel",source="EdwardsLab",type="Subsumption")

    m.set_annotation(property="Subsumes",value=list(subsumes[acc]),source="EdwardsLab",type="Subsumption")
    m.set_annotation(property="SubsumedBy",value=list(subsumedby[acc]),source="EdwardsLab",type="Subsumption")
    
    if w.put(m):
        print acc
