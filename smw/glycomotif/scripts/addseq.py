#!/bin/env python27

from collections import defaultdict

from getwiki import GlycoMotifWiki

import findpygly
from pygly.GlycanResource import GlyTouCan
import sys

w = GlycoMotifWiki()

gtc = None

def accs():
    global gtc
    if len(sys.argv) > 1:
        gtc = GlyTouCan(prefetch=False,usecache=False)
	for acc in sys.argv[1:]:
	    yield acc
    else:
        gtc = GlyTouCan(prefetch=True,usecache=False)
	for m in w.itermotif():
	    yield m.get('id')

for acc in accs():
    m = w.get(acc)
    gtcacc = m.get('glytoucan')
    wurcs = gtc.getseq(gtcacc,'wurcs')
    glycoct = gtc.getseq(gtcacc,'glycoct')
    if not glycoct:
	glycoct = gtc.glycoct(gtcacc)
    if wurcs:
        m.set('wurcs',wurcs)
    if glycoct:
        m.set('glycoct',glycoct)
    if w.put(m):
	print acc
