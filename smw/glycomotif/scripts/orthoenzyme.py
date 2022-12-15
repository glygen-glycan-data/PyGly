#!/bin/env python27

from getwiki import GlycoMotifWiki, Enzyme
import sys, re, glob, json

w = GlycoMotifWiki()

mgi = dict()
ortho = dict()
for l in open(sys.argv[1]):
    sl = l.split()
    if sl[0].upper() == sl[1].upper():
        ortho[sl[0]] = sl[1]
        ortho[sl[1]] = sl[0]
    else:
        if sl[0] not in ortho:
            ortho[sl[0]] = sl[1]
        if sl[1] not in ortho:
            ortho[sl[1]] = sl[0]
    mgi[sl[1]] = sl[2]

for e in w.iterenzyme():
    gn = e.get('genename')
    ogn = ortho.get(gn)
    if not ogn:
        continue
    o = w.get(ogn)
    if o:
        e.set("ortholog",ogn)
    if gn in mgi:
        e.set("mgiacc",mgi[gn])
    if w.put(e):
        print gn
