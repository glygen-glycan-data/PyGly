#!/bin/env python3.12

import sys
from collections import defaultdict

from getwiki import GlycanData, Glycan
w = GlycanData()

import findpygly
from pygly.GlycanResource import MatrixDBSourceFile

mdb = MatrixDBSourceFile()

matrixdb = defaultdict(set)

for mdbid,gtc,dummy,dummy1 in mdb.allgtc():
    matrixdb[gtc].add(mdbid)

for acc in matrixdb.keys():
    m = w.get(acc)
    if not m:
        continue
    if len(matrixdb[acc]) > 0:
        m.set_annotation(value=list(matrixdb[acc]),property="MatrixDB",source="MatrixDB",type="CrossReference")
    else:
        m.delete_annotations(source="MatrixDB",property="MatrixDB",type="CrossReference")
    if w.put(m):
        print(acc)
