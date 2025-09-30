#!/bin/env python3.12

import sys
from operator import itemgetter
from collections import defaultdict

from getwiki import GlycanData, Glycan
w = GlycanData()

from pygly.GlycanResource import GlycoMotif, GlycoMotifDev

from pygly.GNOme import SubsumptionGraph
gnome = SubsumptionGraph()
gnome.loaddata(sys.argv[1])
sys.argv.pop(1)

debug = False
if len(sys.argv) > 1:
    debug = True

def iterglycan():
    if len(sys.argv) > 1:
        for acc in sys.argv[1:]:
            m = w.get(acc)
            if m:
                yield m
    else:
        for m in w.iterglycan():
            yield m

if debug:
    gm = GlycoMotifDev(local=True,prefetch=False,verbose=True,usecache=True)
else:
    gm = GlycoMotifDev(local=True,prefetch=True,verbose=False,usecache=False)

acc2type = defaultdict(set)
for m in iterglycan():
    acc = m.get('accession')
    print("1:",acc,file=sys.stderr)

    byclsid = defaultdict(dict)
    for clsrow in gm.getclass(acc):
        byclsid[clsrow['classid']][clsrow['classlevel']] = clsrow['classname']

    for clsid in byclsid:
        acc2type[acc].add((byclsid[clsid]['Type'],byclsid[clsid].get('Subtype'),"Direct",clsid))

for acc in sorted(acc2type):
    m = w.get(acc)
    print("2:",acc,file=sys.stderr)

    for ann in acc2type[acc]:
        if ann[2] != "Direct":
            continue
        for anc in gnome.ancestors(acc):
            if anc.startswith('G'):
                acc2type[anc].add((ann[0],ann[1],"Subsumption",acc))

for m in iterglycan():
    acc = m.get('accession')
    print("3:",acc,file=sys.stderr)

    m.delete_annotations(type='Classification')

    freqbytype = defaultdict(int)
    for asn in acc2type[acc]:
        if asn[1]:
            freqbytype[asn[0]] += 1
    retain = set()                                                                                                   
    for asn in acc2type[acc]:
        if asn[1]:
            retain.add(asn)
        elif freqbytype[asn[0]] == 0:
            retain.add(asn)
    acc2type[acc] = retain

    subtypes0 = defaultdict(set)
    subtypes1 = defaultdict(set)
    for asn in acc2type[acc]:
        if asn[2] == "Direct":
            subtypes0[asn[:2]].add(asn[3])
        else:
            subtypes1[asn[:2]].add(asn[3])

    subtypecnt = 0
    stbytcnt = defaultdict(int)
    typecnt = 0
    seenst = set()
    seent = set()
    for st in sorted(subtypes0,key=lambda st: min(subtypes0[st])):
        for sid in sorted(subtypes0[st]):
            if st[0] not in seent:
                m.add_annotation(value=st[0], property='GlycanType', 
                                 source='GlycoMotif', type='Classification', 
                                 sourceid=sid)
                typecnt += 1
                seent.add(st[0])
            if st[1]:
                m.add_annotation(value=st[1], property='GlycanSubtype', 
                                 source='GlycoMotif', type='Classification', 
                                 sourceid=sid)
                subtypecnt += 1
                stbytcnt[st[0]] += 1
            seenst.add(st)
            break
    for st in sorted(subtypes1,key=lambda st: min(subtypes1[st])):
        if st in seenst:
            continue
        for sid in sorted(subtypes1[st]):
            if st[0] not in seent:
                m.add_annotation(value=st[0], property='GlycanType', 
                                 source='Subsumption', type='Classification', 
                                 sourceid=sid)
                typecnt += 1
                seent.add(st[0])
            if st[1]:
                m.add_annotation(value=st[1], property='GlycanSubtype', 
                                 source='Subsumption', type='Classification', 
                                 sourceid=sid)
                subtypecnt += 1
                stbytcnt[st[0]] += 1
            break

    m.set_annotation(value=typecnt, property='GlycanTypeCount', source='GlycanData', type='Classification')
    m.set_annotation(value=subtypecnt, property='GlycanSubtypeCount', source='GlycanData', type='Classification')
    m.set_annotation(value=stbytcnt['N-linked'], property='GlycanNLinkedSubtypeCount', source='GlycanData', type='Classification')
    m.set_annotation(value=stbytcnt['O-linked'], property='GlycanOLinkedSubtypeCount', source='GlycanData', type='Classification')
    m.set_annotation(value=stbytcnt['GAG'], property='GlycanGAGSubtypeCount', source='GlycanData', type='Classification')
    m.set_annotation(value=stbytcnt['Glycosphingolipid'], property='GlycanGlycosphingolipidSubtypeCount', source='GlycanData', type='Classification')

    if not debug:
        if w.put(m):
            print("!!",acc,file=sys.stderr)
    else:
        print(m)

