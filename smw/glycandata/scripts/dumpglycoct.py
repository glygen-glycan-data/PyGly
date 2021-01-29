#!/bin/env python27

import sys, time, traceback, hashlib, os, os.path
from collections import defaultdict

import findpygly
from getwiki import GlycanData
w = GlycanData()

for g in w.iterglycan():
    acc = g.get('accession')

    gct = None
    if g.has_annotations(property='GlycoCT',type='Sequence'):
	gct = g.get_annotation_value(property='GlycoCT',type='Sequence')

    comp = None
    if g.has_annotations(property='Composition', type='Structure'):
	comp = g.get_annotation_value(property='Composition', type='Structure')

    undet = None
    if g.has_annotations(property='UndeterminedTopology', type='Structure'):
	undet = g.get_annotation_value(property='UndeterminedTopology', type='Structure')

    nlinked = False
    for an in g.annotations(property='GlycanType',type="Classification"):
	for t in an.get('value'):
	    if t == "N-linked":
		nlinked = True
		break
	if nlinked:
	    break

    nmotif = False
    for ann in g.annotations(property='ClassMotif',type='Motif',source='GlycoMotif'):
        for value in ann.get('value',[]):
            if value == "GGM.001001":
		nmotif = True
		break
	if nmotif:
	    break

    if gct == None or comp!='false' or not nmotif:
	continue

    filename = sys.argv[1] + "/" + acc + ".txt"
    if os.path.exists(filename):
	continue

    print >>sys.stderr, acc
    wh = open(filename,'w')
    wh.write(gct)
    wh.close()

