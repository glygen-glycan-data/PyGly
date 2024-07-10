#!/bin/env python2

import sys
from collections import defaultdict
import csv

from getwiki import GlycanData, Glycan
w = GlycanData()

def pubchemxref(infile):
    d = defaultdict(set)
    f = open(infile,'r')
    # three wacky bytes...
    f.read(3)
    
    for r in csv.DictReader(f):
	cid = r['cid'].strip()
	sid = r['sid'].strip()
	gtcacc = r['sidextid'].strip()
        if gtcacc:
            if cid:
	        d[gtcacc].add("CID"+cid)
            if sid:
	        d[gtcacc].add("SID"+sid)
    return d

pubchem = pubchemxref(sys.argv[1])

for m in w.iterglycan():
    acc = m.get('accession')
    if len(pubchem[acc]) > 0:
        m.set_annotation(value=list(pubchem[acc]), property="PubChem",source="GlyTouCan",type="CrossReference")
    else:
        m.delete_annotations(source="GlyTouCan",property="PubChem",type="CrossReference")
    if w.put(m):
        print acc
