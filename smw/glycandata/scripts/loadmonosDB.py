#!/bin/env python27

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
    monosdic = defaultdict(set)
    acc = g.get('accession')
    glycan = g.getGlycan()
    if not glycan:
        continue
    for m in glycan.all_nodes():
        try:
            glycoctsym = glycoctformat.mtoStr(m)
        except KeyError:
            print acc,"-"
            continue
        try:
            monosdic[acc].add(monosdb[glycoctsym])
        except KeyError:
            continue
    g.delete_annotations(source="EdwardsLab",property="MonosaccharideDB",type="CrossReference")
    g.set_annotation(value=list(monosdic[acc]), property="MonosaccharideDB",source="EdwardsLab",type="CrossReference")
    if w.put(g):
        print acc


