#!/bin/env python27

import sys, time, re
from collections import defaultdict
import csv

gtc2uc = defaultdict(set)
for l in open(sys.argv[1]):
    sl = l.split()
    gtc2uc[sl[1]].add(sl[0])

uc2taxa = defaultdict(set)
for l in open(sys.argv[2]):
    sl = l.split()
    if sl[0].startswith('Hex'):
	sl[0] = 'comp_' + sl[0]
    uc2taxa[sl[0]].add(sl[1])

gtc2name = defaultdict(set)
for fn in sys.argv[3:]:
    for l in open(fn):
	sl = l.split()
	if not re.search(r'^HexNAc\d+Hex\d+dHex\d+NeuAc\d+NeuGc\d+Pent\d+S\d+P\d+KDN\d+HexA\d+$',sl[1]):
	    continue
	compname = "comp_" + sl[1]
	gtc2name[sl[0]].add(compname)

for acc in sorted(set(list(gtc2uc) + list(gtc2name))):
    ucaccs = gtc2uc[acc]
    ucaccs.update(gtc2name[acc])
    for ucacc in ucaccs:
        uctaxa = uc2taxa[ucacc]
	for taxid in uctaxa:
	    print "\t".join(map(str,[acc,taxid,"UniCarbKB",ucacc]))
