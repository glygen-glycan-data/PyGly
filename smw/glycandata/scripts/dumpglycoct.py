#!/bin/env python2

import sys, time, traceback, hashlib, os, os.path, glob
from collections import defaultdict

import findpygly
from getwiki import GlycanData
w = GlycanData()

allfn = set(glob.glob(sys.argv[1]+"/G*.txt"))

glycanclass = sys.argv[2]
assert glycanclass in ("N-linked","O-linked")

for g in w.iterglycan():
    acc = g.get('accession')

    gct = None
    if g.has_annotations(property='GlycoCT',type='Sequence'):
	gct = g.get_annotation_value(property='GlycoCT',type='Sequence')

    inclass = False
    if glycanclass == "N-linked":
        if g.has_annotations(property='ClassMotif',type='Motif',source='GlycoMotif'):
	    for value in g.get_annotation_values(property='ClassMotif',type='Motif',source='GlycoMotif'):
                if value == "GGM.001001":
		    inclass = True
		    break

    elif glycanclass == "O-linked":
        if g.has_annotations(property='ClassMotif',type='Motif',source='GlycoMotif'):
	    for value in g.get_annotation_values(property='ClassMotif',type='Motif',source='GlycoMotif'):
                if value in ("GGM.001006","GGM.001010","GGM.001014","GGM.001016","GGM.001018"):
		    inclass = True
		    break
    else:
	raise RuntimeError("Bad glycan-class...")

    if gct == None or not inclass:
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
