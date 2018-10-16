#!/bin/env python27

import sys, traceback, os, csv

from getwiki import GlycoMotifWiki
from getwiki import GlydinCummingsMotif, GlydinHayesMotif, GlydinCermavMotif, GlydinSugarbindMotif, GlydinBioligoMotif

w = GlycoMotifWiki()

from pygly.GlyTouCan import GlyTouCan
gtc = GlyTouCan()

basedir = sys.argv[1]
classnames = [GlydinCummingsMotif, GlydinHayesMotif, GlydinCermavMotif, GlydinSugarbindMotif, GlydinBioligoMotif]

aglycon2stdaglycon = {
    "Ser/Thr O-Mannose": "Ser/Thr",
    "Ceramide": "Cer",
    "Cer 9-O-Acetyl": "Cer",
}
stdaglycon = ["Ser/Thr","Cer","R","Other"]
reaglycon = ["Ser/Thr","Cer", "Other"]

current = set()
for motifClass in classnames:

    filename = os.path.join(basedir,"Glydin-" +  motifClass.id + ".tsv")
    if not os.path.exists(filename):
	continue
    rows = csv.reader(open(filename),dialect='excel-tab')

    for linenum,r in enumerate(rows):
        if linenum == 0:
            continue

        content = map(lambda x: x.strip(), r)
        
        accession = "%06d" % int(content[0])
        glytoucan = content[1]
        aglycon = content[2]
        redend = False
        name = None
        
        if aglycon:
            if aglycon in aglycon2stdaglycon:
                aglycon = aglycon2stdaglycon[aglycon]
            elif aglycon in stdaglycon:
                pass
            else:
                aglycon = "Other"
            
            if aglycon in reaglycon:
                redend = True
        else:
            aglycon = None

        if len(content) > 3 and content[3]:
            name = content[3]
 
        motif = motifClass(accession=accession, name=name, glytoucan=glytoucan, redend=redend, aglycon=aglycon)

        if w.update(motif):
            print accession
        current.add(accession)


    for m in w.itermotif(collection=motifClass):
        if m.get('accession') not in current:
            print "Deleting:", m.get('pagename')
            w.delete(m.get('pagename'))
