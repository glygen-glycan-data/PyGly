#!/bin/env python27

import sys, time, traceback
from collections import defaultdict
from operator import itemgetter

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
    for acc,gtcacc,alignment,redend,aglycon,names,pmids,keywords in sorted(gm.allmotifs(coll)):
	if coll == "GTC":
	    continue
	if coll == "GGM" and int(acc.split('.')[1]) >= 1000:
	    continue
	if allmotifs:
	    print "\t".join([acc,"Collection",coll])
	    print "\t".join([acc,"Accession",acc.split('.',1)[1]])
	    for i,n in enumerate(names):
		n = n.replace(u'\u2013','-')
	        if i == 0:
	            print "\t".join([acc,"PreferredName",n])
	        else:
	            print "\t".join([acc,"AlternativeName",n])
	    if aglycon:
	        print "\t".join([acc,"Aglycon",aglycon])
	    print "\t".join([acc,"Alignment",alignment])
	    print "\t".join([acc,"MotifGlyTouCan",gtcacc])
	    print "\t".join([acc,"ReducingEnd",redend])
	    for i,pmid in enumerate(pmids):
		print "\t".join([acc,"PMID",pmid])
	    for i,kw in enumerate(keywords):
		print "\t".join([acc,"Keyword",kw])
	allmotifdata[acc] = dict(alignment=alignment,label=names[0])

if allmotifs:
    sys.exit(0)

print "\t".join(["GlyTouCanAccession","MotifAccession","Label","Alignment"])
for g in w.iterglycan():

    acc = g.get('accession')

    motifs = set()
    for coll in collections:
        motifs.update(map(itemgetter(0),filter(lambda t: t[0] in allmotifdata and t[1],gm.getmotif(coll,acc))))

    for motifacc in sorted(motifs):
	print "\t".join([acc,motifacc,allmotifdata[motifacc]['label'],allmotifdata[motifacc]['alignment']])

