#!/bin/env python27

from getwiki import GlycoMotifWiki, GlyTouCanMotif

w = GlycoMotifWiki()

for m in w.itermotif():
    print m.get('collection'), m.get('accession'), m.get('id')
