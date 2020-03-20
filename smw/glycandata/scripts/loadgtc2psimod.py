#!/bin/env python27

import sys
from collections import defaultdict
import csv

from getwiki import GlycanData, Glycan
w = GlycanData()

def xref(infile):
    d = defaultdict(set)
    f = open(infile,'r')
    
    for r in csv.reader(f,dialect='excel-tab'):
	gtcacc = r[0]
	otheracc = r[1]
	d[gtcacc].add(otheracc)
    return d

xrefs = xref(sys.argv[1])

for m in w.iterglycan():
    acc = m.get('accession')
    if len(xrefs[acc]) > 0:
        m.set_annotation(value=list(xrefs[acc]), property="PSI-MOD",source="EdwardsLab",type="CrossReference")
    else:
        m.delete_annotations(source="GlyTouCan",property="PSI-MOD",type="CrossReference")
    if w.put(m):
        print acc
