#!/bin/env python27

import sys,traceback
import findpygly
from pygly.GlycoMotifWiki import GlycoMotifWiki, GlycoEpitopeMotif
w = GlycoMotifWiki()

import csv

current = set()
for r in csv.DictReader(open(sys.argv[1]),dialect='excel-tab'):
    if not r['glytoucan']:
	continue
    motif = GlycoEpitopeMotif(accession=r['acc'],name=r['name'],glytoucan=r['glytoucan'],redend="False")
    if w.put(motif):
	print r['acc']
    current.add(r['acc'])

for m in w.itermotif():
    if m.get('collection') == "GlycoEpitope" and \
       m.get('accession') not in current:
        print "Deleting:",m.get('pagename')
        w.delete(m.get('pagename'))
