#!/bin/env python3.12

import sys
from collections import defaultdict
import csv

from getwiki import GlycanData, Glycan
w = GlycanData()

import findpygly
from pygly.GlycanResource import GlyConnectSourceFile

ggsf = GlyConnectSourceFile()

gtc2acc = defaultdict(set)
for row in ggsf.allgtc():
    if row[1] != "-":
        gtc2acc[row[1]].add(row[0])

for m in w.iterglycan():
    acc = m.get('accession')
    xrefs = set()
    if m.has_annotations(property="GlyConnectStructure", source="GlyTouCan", type="CrossReference"):
        xrefs = set(m.get_annotation_values(property="GlyConnectStructure", source="GlyTouCan", type="CrossReference"))
    m.delete_annotations(property="GlyConnectStructure", source="GlyConnect", type="CrossReference")
    m.delete_annotations(property="GlyConnectComposition", source="GlyConnect", type="CrossReference")
    for gcacc in gtc2acc[acc]:
        if gcacc.startswith("S") and gcacc[1:] not in xrefs:
            m.add_annotation(value=gcacc[1:], property="GlyConnectStructure", source="GlyConnect", type="CrossReference")
        elif gcacc.startswith("C"):
            m.add_annotation(value=gcacc[1:], property="GlyConnectComposition", source="GlyConnect", type="CrossReference")
    if w.put(m):
        print(acc)
