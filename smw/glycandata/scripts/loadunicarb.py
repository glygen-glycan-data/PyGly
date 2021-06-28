#!/bin/env python2

import sys, time, re
from collections import defaultdict
import csv

from getwiki import GlycanData, Glycan
w = GlycanData()

gtc2uc = defaultdict(set)
for l in open(sys.argv[1]):
    sl = l.split()
    gtc2uc[sl[1]].add(sl[0])

uc2pubmed = defaultdict(set)
for l in open(sys.argv[2]):
    sl = l.split()
    if sl[0].startswith('Hex'):
	sl[0] = 'comp_' + sl[0]
    uc2pubmed[sl[0]].add(sl[1])

gtc2name = defaultdict(set)
for fn in sys.argv[3:]:
    for l in open(fn):
	sl = l.split()
	if not re.search(r'^HexNAc\d+Hex\d+dHex\d+NeuAc\d+NeuGc\d+Pent\d+S\d+P\d+KDN\d+HexA\d+$',sl[1]):
	    continue
	compname = "comp_" + sl[1]
	gtc2name[sl[0]].add(compname)

for m in w.iterglycan():
    start = time.time()
    acc = m.get('accession')

    ucaccs = gtc2uc[acc]
    ucaccs.update(gtc2name[acc])
    # Remove any GlyTouCan xref to UniCarbKB...
    m.delete_annotations(property="UniCarbKB",type="CrossReference")
    if len(ucaccs) > 0:
	m.set_annotation(value=list(ucaccs), property="UniCarbKB", source="UniCarbKB", type="CrossReference")

    # m.delete_annotations(property="Taxonomy", source="UniCarbKB", type="Taxonomy")
    # for ucacc in ucaccs:
    #	uctaxa = uc2taxa[ucacc]
    #     if len(uctaxa) > 0:
    # 	    m.set_annotation(value=list(uctaxa), property="Taxonomy", source="UniCarbKB", sourceid=ucacc, type="Taxonomy")

    m.delete_annotations(property="Publication", source="UniCarbKB", type="Publication")
    for ucacc in ucaccs:
	ucpubmed = uc2pubmed[ucacc]
        if len(ucpubmed) > 0:
	    m.set_annotation(value=list(ucpubmed), property="Publication", source="UniCarbKB", sourceid=ucacc, type="Publication")

    if w.put(m):
        print >>sys.stderr, "%s updated in %.2f sec"%(m.get('accession'),time.time()-start,)
    else:
        print >>sys.stderr, "%s checked in %.2f sec"%(m.get('accession'),time.time()-start,)
