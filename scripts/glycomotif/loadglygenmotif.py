#!/bin/env python27

from collections import defaultdict

import findpygly

from pygly.GlycoMotifWiki import GlycoMotifWiki, GlyGenMotif
w = GlycoMotifWiki()

from pygly.GlyTouCan import GlyTouCan
gtc = GlyTouCan()

gtc2motif = defaultdict(list)
for m in w.itermotif():
    if m.get('collection') in ('UGA-CCRC','GlyTouCan','GlycoEpitope'):
        gtc = m.get('glytoucan')
        gtc2motif[gtc].append(m)

for gtc,motiflist in sorted(gtc2motif.iteritems()):
    names = []
    redend = set()
    for m in motiflist:
	if m.has('name'):
	    names.extend(m.get('name'))
	if m.has('redend'):
	    redend.add(m.get('redend',False))
    if len(redend) == 1:
	redend = redend.pop()
    else:
	redend = False
    motif = GlyGenMotif(accession=gtc,name=names,redend=redend)
    if w.update(motif):
	print gtc
