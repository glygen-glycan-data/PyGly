#!/bin/env python27

import sys,traceback

from getwiki import GlycoMotifWiki, CCRCMotif
w = GlycoMotifWiki()

from pygly.GlyTouCan import GlyTouCan
gtc = GlyTouCan()

from gtccache import GlyTouCanCache
gtccache = GlyTouCanCache()

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

    redend = None
    aglycon = None
    redendstr = r.get('end')
    if redendstr != None:
	redendstr = redendstr.strip()
    if redendstr:
        if redendstr in ("-Cer","Cer"):
	    redend = True
	    aglycon = "Cer"
        elif redendstr in ("-Ser/Thr","Ser/Thr"):
	    redend = True
	    aglycon = "Ser/Thr"
        elif redendstr in ("-R","R"):
	    redend = False
	    aglycon = "R"
	
    if not glycoct and not glytoucan:
	continue

    accession = "%06d"%(index,)

    if not glytoucan:
	glytoucan = gtccache.id2gtc(CCRCMotif.id + "." + accession)
    glycoct1 = gtccache.gtc2glycoct(glytoucan)
    wurcs = gtccache.gtc2wurcs(glytoucan)

    if not glytoucan:
	try:
	    glytoucan,isnew = gtc.register(glycoct)
	except:
	    # traceback.print_exc()
	    continue

    motif = CCRCMotif(accession=accession,name=name,glytoucan=glytoucan,redend=redend,aglycon=aglycon,glycoct=glycoct1,wurcs=wurcs)
    if not motif.get('glycoct'):
        motif.set('glycoct',glycoct)
    if w.put(motif):
	print >>sys.stderr, accession
    current.add(accession)

for m in w.itermotif(collection=CCRCMotif):
    if m.get('accession') not in current:
        print >>sys.stderr, "Deleting:",m.get('pagename')
        w.delete(m.get('pagename'))
