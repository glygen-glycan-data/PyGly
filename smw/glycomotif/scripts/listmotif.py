#!/bin/env python27

from getwiki import GlycoMotifWiki, GlyTouCanMotif

w = GlycoMotifWiki()
print w

for m in w.itercat('Motif'):
    print m
    break

for m in w.itermotif():
    print m
    break
