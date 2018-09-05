#!/bin/env python27

import findpygly
from pygly.GlycoMotifWiki import GlycoMotifWiki, GlyTouCanMotif

w = GlycoMotifWiki()

from pygly.GlyTouCan import GlyTouCan

gtc = GlyTouCan()

current = set()
for m,l,re in sorted(gtc.allmotifs()):
    redend = ("True" if re=="true" else "False")
    motif = GlyTouCanMotif(accession=m,name=l,redend=redend)
    if w.put(motif):
	print m
    current.add(m)

for m in w.itermotif():
    if m.get('collection') == "GlyTouCan" and \
       m.get('accession') not in current:
        print "Deleting:",m.get('pagename')
        w.delete(m.get('pagename'))
