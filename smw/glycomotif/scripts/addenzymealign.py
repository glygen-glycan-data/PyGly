#!/bin/env python27

from getwiki import GlycoMotifWiki
import sys, re, csv, gzip
from collections import defaultdict

w = GlycoMotifWiki()

alignmap = dict(core="Core",
                substr="Substructure",
                nonred="Nonreducing-End",
                whole="Whole-Glycan")

enzdata = defaultdict(set)
laststruct = None
for r in csv.DictReader(gzip.open(sys.argv[1]),dialect='excel-tab'):
    mgtc = r['Motif']
    if r['MotifResidue'] == '-':
        continue
    if r['Structure'] != laststruct:
        print(r['Structure'])
        laststruct = r['Structure']
    align = alignmap[r['AlignmentType']]
    he = r['HumanEnzyme'].split(',')
    me = r['MouseEnzyme'].split(',')
    enz = set(he + me)
    enz = list(filter(lambda e: e != "-",enz))
    enzdata[mgtc,align].update(enz)

for m in w.itermotif(collection="GGM"):
    mgtc = m.get('glytoucan')
    align = m.get('alignment')[0]
    enzymes = m.get("enzyme",set())
    if '-' in enzymes:
        enzymes.remove('-')
    enzymes.update(enzdata[mgtc,align])
    if len(enzymes) > 0:
        m.set("enzyme",enzymes)
    else:
        m.delete("enzyme")
    if w.put(m):
        print >>sys.stderr, m.get('id')
