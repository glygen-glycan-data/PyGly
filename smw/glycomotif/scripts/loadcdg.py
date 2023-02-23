#!/bin/env python2

from getwiki import GlycoMotifWiki, Enzyme
import sys, re, glob, json

w = GlycoMotifWiki()

cdg = set(open(sys.argv[1]).read().split())

for e in w.iterenzyme():
    gn = e.get('genename')
    if gn in cdg:
       e.set("iscdg",'true')
    else:
       e.delete("iscdg")
    if w.put(e):
        print gn
