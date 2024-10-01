#!/bin/env python3.12

import sys, time, traceback
from collections import defaultdict
from operator import itemgetter

from getwiki import GlycanData
w = GlycanData()

import findpygly
from pygly.GlycanResource import GlycoMotif

gm = GlycoMotif(prefetch=True,verbose=False,usecache=False)

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
    print("\t".join(["MotifAccession","Property","Value"]))
allmotifdata = dict()
for coll in collections:
    for acc,gtcacc,alignment,redend,aglycon,names,pmids,keywords,dbxrefs in sorted(gm.allmotifs(coll)):
        if coll == "GTC":
            continue
        if coll == "GGM" and int(acc.split('.')[1]) >= 1000:
            continue
        if allmotifs:
            print("\t".join([acc,"Collection",coll]))
            print("\t".join([acc,"Accession",acc.split('.',1)[1]]))
            names = [  n.replace(u'\u2013','-') for n in names ]
            print("\t".join([acc,"PreferredName",names[0]]))
            for n in sorted(names[1:]):
                print("\t".join([acc,"AlternativeName",n]))
            if aglycon:
                print("\t".join([acc,"Aglycon",aglycon]))
            print("\t".join([acc,"Alignment",alignment]))
            print("\t".join([acc,"MotifGlyTouCan",gtcacc]))
            print("\t".join([acc,"ReducingEnd",redend]))
            for i,pmid in enumerate(sorted(pmids,key=int)):
                print("\t".join([acc,"PMID",pmid]))
            for i,kw in enumerate(sorted(keywords)):
                print("\t".join([acc,"Keyword",kw]))
            for i,xr in enumerate(sorted(dbxrefs)):
                if 'GlycanDictionary:' in xr:
                    xr = xr.replace("GlycanDictionary:GD0","GlycanDictionary:GSD0")
                print("\t".join([acc,"CrossRef",xr]))
        allmotifdata[acc] = dict(alignment=alignment,label=names[0])

if allmotifs:
    sys.exit(0)

print("\t".join(["GlyTouCanAccession","MotifAccession","Label","Alignment"]))
for g in w.iterglycan():

    acc = g.get('accession')

    motifs = set()
    for coll in collections:
        motifs.update(map(itemgetter(0),filter(lambda t: t[0] in allmotifdata and t[1],gm.getmotif(coll,acc))))

    for motifacc in sorted(motifs):
        print("\t".join([acc,motifacc,allmotifdata[motifacc]['label'],allmotifdata[motifacc]['alignment']]))

