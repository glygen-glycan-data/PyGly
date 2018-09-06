#!/bin/env python27

from collections import defaultdict

import findpygly

from pygly.GlycoMotifWiki import GlycoMotifWiki, GlyGenMotif
w = GlycoMotifWiki()

from pygly.GlyTouCan import GlyTouCan
gtc = GlyTouCan()

gtc2motif = defaultdict(list)
for m in w.itermotif():
    if m.get('collection') in ('UGA-CCRC','GlyTouCan','GlyGen','GlycoEpitope'):
        gtc = m.get('glytoucan')
        gtc2motif[gtc].append(m)

for gtc,motiflist in sorted(gtc2motif.iteritems()):
    for m in motiflist:
	sameas = map(lambda m1: m1.get('pagename'),motiflist)
	sameas.remove(m.get('pagename'))
	m.set('sameas',sameas)
	w.put(m)
