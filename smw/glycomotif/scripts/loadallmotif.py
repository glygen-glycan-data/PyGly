#!/bin/env python27

from collections import defaultdict

from getwiki import GlycoMotifWiki, AllMotif
w = GlycoMotifWiki()

gtc2motif = defaultdict(list)
for m in w.itermotif():
    if m.get('collection') != AllMotif.id:
        gtc = m.get('glytoucan')
        gtc2motif[gtc].append(m)

for gtc,motiflist in sorted(gtc2motif.iteritems()):
    names = []
    redend = set()
    aglycon = set()
    for m in motiflist:
	if m.has('name'):
	    for n in m.get('name'):
		if n not in names:
	            names.append(n)
	if m.has('redend'):
	    redend.add(m.get('redend',False))
	if m.has('aglycon'):
	    aglycon.update(m.get('aglycon',[]))
    if len(redend) == 1:
	redend = redend.pop()
    else:
	redend = False
    if len(aglycon) > 0:
	aglycon = sorted(aglycon)
    else:
	aglycon = None
    motif = AllMotif(accession=gtc,name=names,redend=redend,aglycon=aglycon)
    if w.update(motif):
	print gtc
