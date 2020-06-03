#!/bin/env python27

import sys, time, traceback
from collections import defaultdict

from getwiki import GlycanData
w = GlycanData()

import findpygly
from pygly.GlycanResource import GlycoMotif

gm = GlycoMotif(prefetch=True,verbose=False)

if sys.argv[1] == "allmotifs":
    allmotifs = True
    allaligns = False
elif sys.argv[1] == "allaligns":
    allmotifs = False
    allaligns = True
else:
    raise RuntimeError("Bad report selection.")

collections = ('GGM','GTC')

if allmotifs:
    print "\t".join(["MotifAccession","Property","Value"])
allmotifdata = dict()
for coll in collections:
    for acc,gtcacc,redend,aglycon,names in sorted(gm.allmotifs(coll)):
	if coll == "GTC":
	    continue
	if coll == "GGM" and int(acc.split('.')[1]) >= 1000:
	    continue
	if allmotifs:
	    print "\t".join([acc,"Collection",coll])
	    print "\t".join([acc,"Accession",acc.split('.',1)[1]])
	    for i,n in enumerate(names):
	        if i == 0:
	            print "\t".join([acc,"PreferredName",n])
	        else:
	            print "\t".join([acc,"AlternativeName",n])
	    if aglycon:
	        print "\t".join([acc,"Aglycon",aglycon])
	    print "\t".join([acc,"MotifGlyTouCan",gtcacc])
	    print "\t".join([acc,"ReducingEnd",redend])
	allmotifdata[acc] = dict(redend=redend,label=names[0])

if allmotifs:
    sys.exit(0)

print "\t".join(["GlyTouCanAccession","MotifAccession","Label","IsReducingEnd"])
for g in w.iterglycan():

    acc = g.get('accession')

    motifs = set()
    for coll in collections:
        motifs.update(filter(lambda a: a in allmotifdata,gm.getmotif(coll,acc)))

    for motifacc in sorted(motifs):
	print "\t".join([acc,motifacc,allmotifdata[motifacc]['label'],allmotifdata[motifacc]['redend']])

