#!/bin/env python27

import sys,traceback

from getwiki import GlycoMotifWiki, GlycoEpitopeMotif
w = GlycoMotifWiki()

import csv
from gtccache import GlyTouCanCache
gtccache = GlyTouCanCache()

current = set()
for r in csv.DictReader(open(sys.argv[1]),dialect='excel-tab'):
    if not r['glytoucan']:
	continue
    redend = None
    aglycon = None
    if r['sequence'].endswith('Cer'):
	redend = True
	aglycon = 'Cer'
    elif r['sequence'].endswith('Ser/Thr'):
	redend = True
	aglycon = 'Ser/Thr'
    elif r['sequence'].endswith('-R'):
	redend = False
	aglycon = 'R'
    wurcs=gtccache.gtc2wurcs(r['glytoucan'])
    glycoct=gtccache.gtc2glycoct(r['glytoucan'])
    motif = GlycoEpitopeMotif(accession=r['acc'],name=r['name'],glytoucan=r['glytoucan'],redend=redend,aglycon=aglycon,wurcs=wurcs,glycoct=glycoct)
    if w.put(motif):
	print >>sys.stderr, r['acc']
    current.add(r['acc'])

for m in w.itermotif(collection=GlycoEpitopeMotif):
    if m.get('accession') not in current:
        print >>sys.stderr, "Deleting:",m.get('pagename')
        w.delete(m.get('pagename'))
