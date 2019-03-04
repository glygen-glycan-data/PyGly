#!/bin/env python27

import sys

from getwiki import GlycanDataWiki, Glycan
w = GlycanDataWiki()

import findpygly
from pygly.GlyTouCan import GlyTouCan

gtc = GlyTouCan()

current = set()
for gtcacc in open(sys.argv[1]):
    gtcacc = gtcacc.strip()
    g = Glycan(accession=gtcacc,mw=gtc.getmass(gtcacc),iupac=gtc.getseq(gtcacc,'iupac_extended'))
    if w.update(g):
	print >>sys.stderr, g.get('accession')
    current.add(gtcacc)

for m in w.iterglycan():
    if m.get('accession') not in current:
        print >>sys.stderr, "Deleting:",m.get('pagename')
        w.delete(m.get('pagename'))
