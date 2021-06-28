#!/bin/env python2

import sys, time, traceback, hashlib, os, os.path, glob
from collections import defaultdict

import findpygly
from getwiki import GlycanData
w = GlycanData()

allfn = set(glob.glob(sys.argv[1]+"/G*.txt"))

for g in w.iterglycan():
    acc = g.get('accession')

    gct = None
    if g.has_annotations(property='GlycoCT',type='Sequence'):
	gct = g.get_annotation_value(property='GlycoCT',type='Sequence')

    nmotif = False
    if g.has_annotations(property='ClassMotif',type='Motif',source='GlycoMotif'):
	for value in g.get_annotation_values(property='ClassMotif',type='Motif',source='GlycoMotif'):
            if value == "GGM.001001":
		nmotif = True
		break

    if gct == None or not nmotif:
	continue

    filename = sys.argv[1] + "/" + acc + ".txt"

    if filename in allfn:
        allfn.remove(filename)	
	continue

    print >>sys.stderr, acc
    wh = open(filename,'w')
    wh.write(gct)
    wh.close()

for fn in allfn:
    print >>sys.stderr, "Removing:",acc
    os.unlink(fn)
