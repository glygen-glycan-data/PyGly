#!/bin/env python2

import sys
from collections import defaultdict

from getwiki import GlycanData, Glycan
from pygly.GlycanFormatter import GlycoCTFormat

w = GlycanData()
glycoctformat = GlycoCTFormat()

monosdb = {}
f = open(sys.argv[1],'r')
for line in f:
    k,v = line.split()
    monosdb[k] = v

for g in w.iterglycan():
    acc = g.get('accession')
    monodbids = set()
    glycan = g.getGlycan()
    if not glycan:
        continue
    for m in glycan.all_nodes():
        try:
            glycoctsym = glycoctformat.mtoStr(m)
        except KeyError:
            continue
        try:
            monodbids.add(monosdb[glycoctsym])
        except KeyError:
            continue
    g.delete_annotations(source="EdwardsLab",property="MonosaccharideDB",type="CrossReference")
    g.set_annotation(value=monodbids, property="MonosaccharideDB",source="EdwardsLab",type="CrossReference")
    if w.put(g):
        print acc

