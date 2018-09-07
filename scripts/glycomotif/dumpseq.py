#!/bin/env python27

import findpygly
from pygly.GlycoMotifWiki import GlycoMotifWiki

import sys, os, os.path

dir = sys.argv[1]
try:
    os.makedirs(dir)
except OSError:
    pass
format = sys.argv[2]
assert format in ('wurcs','glycoct')

seen=set()
w = GlycoMotifWiki()
for m in w.itermotif():
    glytoucan = m.get('glytoucan')
    if glytoucan in seen:
	continue
    seen.add(glytoucan)
    seq = m.get(format)
    if not seq:
	continue
    print glytoucan
    fn = os.path.join(dir,glytoucan+'.txt')
    wh = open(fn,'w')
    wh.write(seq.strip()+"\n")
    wh.close()
