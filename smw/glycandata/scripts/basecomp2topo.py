#!/bin/env python27

import sys

import findpygly

from pygly.GNOme import GNOme
from pygly.GlyTouCan import GlyTouCan

g = GNOme()
gtc = GlyTouCan(usecache=True)

for acc in sys.argv[1:]:
    taccs = set(filter(g.istopology,g.descendants(acc)))
    # print " ".join(sorted(taccs))
    toremove = set()
    for tacc in taccs:
	gly = gtc.getGlycan(tacc)
	if gly.undetermined():
	    toremove.add(tacc)
	for m in gly.all_nodes():
	    if m.stem() == None:
		toremove.add(tacc)
		break
    for tacc in toremove:
        taccs.remove(tacc)
    toremove = set()
    for tacc in taccs:
	if (g.ancestors(tacc) & taccs):
	    toremove.add(tacc)
    for tacc in toremove:
        taccs.remove(tacc)
    print "%s:"%(acc,),",".join(sorted(taccs))
