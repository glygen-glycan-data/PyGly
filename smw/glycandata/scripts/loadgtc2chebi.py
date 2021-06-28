#!/bin/env python2

import sys
from collections import defaultdict
import csv

from getwiki import GlycanData, Glycan
w = GlycanData()

def chebixref(infile):
    d = defaultdict(set)
    f = open(infile,'r')
    
    for r in csv.DictReader(f,dialect='excel-tab'):
	gtcacc = r['gtc'].rsplit('/',1)[-1]
	chebi = r['ch'].rsplit('=',1)[-1]
	d[gtcacc].add(chebi)
    return d

chebi = chebixref(sys.argv[1])

for m in w.iterglycan():
    acc = m.get('accession')
    if len(chebi[acc]) > 0:
        m.set_annotation(value=list(chebi[acc]), property="ChEBI",source="GlyTouCan",type="CrossReference")
    else:
        m.delete_annotations(source="GlyTouCan",property="ChEBI",type="CrossReference")
    if w.put(m):
        print acc
