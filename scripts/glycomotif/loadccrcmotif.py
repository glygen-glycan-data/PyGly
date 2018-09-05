#!/bin/env python27

import sys,traceback
import findpygly
from pygly.GlycoMotifWiki import GlycoMotifWiki, CCRCMotif
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

    redend = r['end']
    if redend:
	redend = redend.strip()
	redend = ("False" if redend == "-R" else "True")
	
    if not glycoct and not glytoucan:
	continue

    if not glytoucan:
	try:
	    glytoucan,isnew = gtc.register(glycoct)
	except:
	    # traceback.print_exc()
	    continue

    accession = "%06d"%(index,)
    motif = CCRCMotif(accession=accession,name=name,glytoucan=glytoucan,redend=redend)
    if w.put(motif):
	print "%06d"%(index,)
    current.add(accession)

for m in w.itermotif():
    if m.get('collection') == "CCRC" and \
       m.get('accession') not in current:
        print "Deleting:",m.get('pagename')
        w.delete(m.get('pagename'))
