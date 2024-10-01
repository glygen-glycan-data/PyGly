#!/bin/env python3.12

import sys
from collections import defaultdict
import csv

from getwiki import GlycanData, Glycan
w = GlycanData()

def mdbxref(infile):
    d = defaultdict(set)
    f = open(infile,'r')
    
    for r in csv.DictReader(f,dialect='excel-tab'):
        gtcacc = None
        for exid in r['External identifiers'].split(','):
            if exid.startswith('GlyTouCan:'):
                gtcacc = exid.split(':',1)[1]
                break
        if not gtcacc:
            continue
        mdbid = r['MatrixDB identifier'].strip()
        if not mdbid.startswith('GAG_'):
            continue
        d[gtcacc].add(mdbid)
    return d

matrixdb = mdbxref(sys.argv[1])

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
