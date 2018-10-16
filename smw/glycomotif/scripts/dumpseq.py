#!/bin/env python27

import findpygly
from getwiki import GlycoMotifWiki
# we must instantiate early, before any use of the command-line. 
w = GlycoMotifWiki()
#print w._prefix

import sys, os, os.path

dir = sys.argv[1]
try:
    os.makedirs(dir)
except OSError:
    pass
format = sys.argv[2]
assert format in ('wurcs','glycoct')

seen=set()
for m in w.itermotif():
    glytoucan = m.get('glytoucan')
    if glytoucan in seen:
        continue
    seen.add(glytoucan)
    seq = m.get(format)
    if not seq:
        continue
    print glytoucan
    # fn = os.path.join(dir,glytoucan+'.txt')
    # wh = open(fn,'w')
    # wh.write(seq.strip()+"\n")
    # wh.close()
