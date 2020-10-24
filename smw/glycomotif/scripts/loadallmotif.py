#!/bin/env python27

from collections import defaultdict
import sys

from getwiki import GlycoMotifWiki, AllMotif
w = GlycoMotifWiki()

import findpygly
from pygly.GlycanResource import GlyTouCan

gtc2motif = defaultdict(list)
for m in w.itermotif():
    print(m.get('id'))
    if m.get('collection') != AllMotif.id:
        gtc = m.get('glytoucan')
        gtc2motif[gtc].append(m)
print('Done')

current = set()
for gtc,motiflist in sorted(gtc2motif.iteritems()):
    print gtc
    names = []
    redend = set()
    aglycon = set()
    alignment = set()
    for m in sorted(motiflist,key=lambda m: (m.get('collection'),m.get('accession'))):
	if m.has('name'):
	    for n in m.get('name'):
		if n not in names:
	            names.append(n)
	if m.has('redend'):
	    redend.update(m.get('redend',[]))
	if m.has('aglycon'):
	    aglycon.update(m.get('aglycon',[]))
	if m.has('alignment'):
	    alignment.update(m.get('alignment',[]))
    if len(redend) == 0:
	redend = None
    if len(aglycon) == 0:
	aglycon = None
    if len(alignment) == 0:
	alignment = None
    motif = AllMotif(accession=gtc,name=names,redend=redend,aglycon=aglycon,alignment=alignment)
    current.add(gtc)
    if w.update(motif):
	print >>sys.stderr, gtc

for m in w.itermotif(collection=AllMotif):
    if m.get('accession') not in current:
        print >>sys.stderr, "Deleting:",m.get('pagename')
        w.delete(m.get('pagename'))
