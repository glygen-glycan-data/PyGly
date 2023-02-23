#!/bin/env python27

from getwiki import GlycoMotifWiki, Enzyme
import sys, re, glob, json
from collections import defaultdict

w = GlycoMotifWiki()

mgi = dict()
ortho = defaultdict(set)
for l in open(sys.argv[1]):
    sl = l.split()
    ortho[sl[0]].add(sl[1])
    ortho[sl[1]].add(sl[0])
    mgi[sl[1]] = sl[2]

for e in w.iterenzyme():
    gn = e.get('genename')
    goodogn = set()
    for ogn in ortho.get(gn,[]):
      o = w.get(ogn)
      if o:
          goodogn.add(ogn)
    if len(goodogn) > 0:
        e.set("ortholog",goodogn)
    else:
        e.delete("ortholog")
    if gn in mgi:
        e.set("mgiacc",mgi[gn])
    if w.put(e):
        print gn
