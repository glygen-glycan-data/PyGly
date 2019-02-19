#!/bin/env python27

from collections import defaultdict

from getwiki import GlycoMotifWiki

w = GlycoMotifWiki()

gtc2motif = defaultdict(list)
for m in w.itermotif():
    gtc = m.get('glytoucan')
    if gtc:
        gtc2motif[gtc].append(m)

for gtc,motiflist in sorted(gtc2motif.iteritems()):
    print gtc
    for m in motiflist:
	sameas = sorted(map(lambda m1: m1.get('id'),motiflist))
	sameas.remove(m.get('id'))
	m.set('sameas',sameas)
	w.update(m)
