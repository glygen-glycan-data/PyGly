#!/bin/env python27

import sys

from getwiki import GlycoMotifWiki, GlyTouCanMotif
w = GlycoMotifWiki()

import findpygly
from pygly.GlyTouCan import GlyTouCan

gtc = GlyTouCan()

from gtccache import GlyTouCanCache
gtccache = GlyTouCanCache()

current = set()
for m,l,re in sorted(gtc.allmotifs()):
    motif = GlyTouCanMotif(accession=m,name=l,redend=re,
			   wurcs = gtccache.gtc2wurcs(m),
			   glycoct = gtccache.gtc2glycoct(m))
    if w.put(motif):
	print >>sys.stderr, m
    current.add(m)

for m in w.itermotif(collection=GlyTouCanMotif):
    if m.get('accession') not in current:
        print >>sys.stderr, "Deleting:",m.get('pagename')
        w.delete(m.get('pagename'))
