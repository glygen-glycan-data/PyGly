#!/bin/env python27

from getwiki import GlycoMotifWiki, GlyTouCanMotif
w = GlycoMotifWiki()

import findpygly
from pygly.GlyTouCan import GlyTouCan

gtc = GlyTouCan()

current = set()
for m,l,re in sorted(gtc.allmotifs()):
    motif = GlyTouCanMotif(accession=m,name=l,redend=re)
    if w.update(motif):
	print m
    current.add(m)

for m in w.itermotif(collection=GlyTouCanMotif):
    if m.get('accession') not in current:
        print >>sys.stderr, "Deleting:",m.get('pagename')
        w.delete(m.get('pagename'))
