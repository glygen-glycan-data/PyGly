#!/bin/env python27

import sys,traceback

from getwiki import GlycoMotifWiki, GlyGenMotif
w = GlycoMotifWiki()

from dataset import TSVFileTable
rows = TSVFileTable(sys.argv[1])

current = set()
for r in rows:

    accession = "%06d"%(int(r['accession']),)
    name = [r['name'].strip()]
    gtc = r['gtc'].strip()
    redend = r['redend'].strip()
    aglycon = r['aglycon'].strip()
    for h in r.keys():
	if h.startswith('altname'):
	    name.append(r[h].strip())

    motif = GlyGenMotif(accession=accession,name=name,glytoucan=gtc,redend=redend,aglycon=aglycon)
    if w.update(motif):
	print >>sys.stderr, accession
    current.add(accession)

for m in w.itermotif(collection=GlyGenMotif):
    if m.get('accession') not in current:
        print >>sys.stderr, "Deleting:",m.get('pagename')
        w.delete(m.get('pagename'))
