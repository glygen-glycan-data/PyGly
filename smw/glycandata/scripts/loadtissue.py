#!/bin/env python3.12

import sys
from operator import itemgetter
from collections import defaultdict

from getwiki import GlycanData, Glycan
w = GlycanData()

from pygly.GlycanResource import GlyGenSourceFile

from pygly.GNOme import SubsumptionGraph
gnome = SubsumptionGraph()
gnome.loaddata(sys.argv[1])
sys.argv.pop(1)

ggsf = GlyGenSourceFile()

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

acc2tissue = defaultdict(set)
for r in ggsf.allsourcetissue():
    d = dict(zip("sid,gtc,taxid,tissue,source,sourceid".split(','),r))
    acc = d.get('gtc')
    print("1:",acc,file=sys.stderr)
    if d['source'] == 'GlyCosmos':
        d['sourceid'] = acc
    acc2tissue[acc].add((d['taxid'],d['tissue'],"Direct",d['source'],d['sourceid']))


for acc in sorted(acc2tissue):
    print("2:",acc,file=sys.stderr)

    for ann in acc2tissue[acc]:
        if ann[2] != "Direct":
            continue
        for anc in gnome.ancestors(acc):
            if anc.startswith('G'):
                acc2tissue[anc].add((ann[0],ann[1],"Subsumption","GNOme",acc))

taxid2prop = {
    '9606': 'TissueInHuman',
    '7227': 'TissueInFruitFly',
    '10090': 'TissueInMouse',
    '7955': 'TissueInZebraFish',
    '9913': 'TissueInBovine',
    '9823': 'TissueInPig',
    '9031': 'TissueInChicken',
    '10116': 'TissueInRat',
}

for m in iterglycan():
    acc = m.get('accession')
    print("3:",acc,file=sys.stderr)

    m.delete_annotations(type='Tissue')

    tissues = defaultdict(list)
    for t in sorted(acc2tissue[acc],key=lambda t: (t[2] == "Direct")*(-1)):
        if t[:2] not in tissues:
            tissues[t[:2]].append(t[2:])
        else:
            hasdirect = tissues[t[:2]][0][0] == "Direct"
            if hasdirect:
                 if t[2] == "Direct":
                     tissues[t[:2]].append(t[2:])
            else:
                tissues[t[:2]].append(t[2:])

    for k,vl in tissues.items():
        taxid,tis=k
        prop = taxid2prop[taxid]
        for v in vl:
            m.add_annotation(value=tis,property=prop,
                             source=v[1],type='Tissue',sourceid=v[2])

    if not debug:
        if w.put(m):
            print("!!",acc,file=sys.stderr)
    else:
        print(m)

