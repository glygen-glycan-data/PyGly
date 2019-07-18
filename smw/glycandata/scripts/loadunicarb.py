#!/bin/env python27

import sys, time
from collections import defaultdict
import csv

from getwiki import GlycanData, Glycan
w = GlycanData()

gtc2uc = defaultdict(set)
for l in open(sys.argv[1]):
    sl = l.split()
    gtc2uc[sl[1]].add(sl[0])

uc2taxa = defaultdict(set)
for taxafile in sys.argv[2:]:
  for l in open(taxafile):
    sl = l.split()
    uc2taxa[sl[0]].add(sl[1])

for m in w.iterglycan():
    start = time.time()
    acc = m.get('accession')

    ucaccs = gtc2uc[acc]
    # Remove any GlyTouCan xref to UniCarbKB...
    m.delete_annotations(property="UniCarbKB",type="CrossReference")
    if len(ucaccs) > 0:
	m.set_annotation(value=list(ucaccs), property="UniCarbKB", source="UniCarbKB", type="CrossReference")

    uctaxa = set()
    for ucacc in ucaccs:
	uctaxa.update(uc2taxa[ucacc])
    if len(uctaxa) > 0:
	m.set_annotation(value=list(uctaxa), property="Taxonomy", source="UniCarbKB", type="Taxonomy")
    else:
        m.delete_annotations(property="Taxonomy", source="UniCarbKB", type="Taxonomy")

    if w.put(m):
        print >>sys.stderr, "%s updated in %.2f sec"%(m.get('accession'),time.time()-start,)
    else:
        print >>sys.stderr, "%s checked in %.2f sec"%(m.get('accession'),time.time()-start,)
