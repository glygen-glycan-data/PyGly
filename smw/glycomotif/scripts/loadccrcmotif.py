#!/bin/env python27

import sys,traceback

from getwiki import GlycoMotifWiki, CCRCMotif
w = GlycoMotifWiki()

from pygly.GlyTouCan import GlyTouCan
gtc = GlyTouCan()

from dataset import XLSXFileTable
rows = XLSXFileTable(sys.argv[1])

current = set()
for r in rows:

    try:
        index = int(r['#'])
    except ValueError:
	traceback.print_exc()
	continue

    name = r['Trivial Name-Cummings']
    if name:
	name = name.strip()

    if not name:
	continue

    glycoct = r['Glyco CT']
    if glycoct:
	glycoct = glycoct.strip()

    glytoucan = r['GlyToucan ID']
    if glytoucan:
	glytoucan = glytoucan.strip()

    redend = False
    aglycon = None
    redendstr = r['end']
    if redendstr != None:
	redendstr = redendstr.strip()
        if redendstr in ("-Cer","Cer"):
	    redend = True
	    aglycon = "Cer"
        elif redendstr in ("-Ser/Thr","Ser/Thr"):
	    redend = True
	    aglycon = "Ser/Thr"
        elif redendstr in ("-R","R"):
	    aglycon = "R"
	
    if not glycoct and not glytoucan:
	continue

    if not glytoucan:
	try:
	    glytoucan,isnew = gtc.register(glycoct)
	except:
	    # traceback.print_exc()
	    continue

    accession = "%06d"%(index,)
    motif = CCRCMotif(accession=accession,name=name,glytoucan=glytoucan,redend=redend,aglycon=aglycon)
    if w.update(motif):
	print accession
    current.add(accession)

for m in w.itermotif(collection=CCRCMotif):
    if m.get('accession') not in current:
        print "Deleting:",m.get('pagename')
        w.delete(m.get('pagename'))
