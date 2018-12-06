#!/bin/env python27

import findpygly
from getwiki import GlycoMotifWiki, AllMotif
# we must instantiate early, before any use of the command-line.
w = GlycoMotifWiki()


import sys, os, os.path

dir = "dumps"
try:
    os.makedirs(dir)
except OSError:
    pass

try:
    os.makedirs(os.path.join(dir, "glycoct"))
except OSError:
    pass

try:
    os.makedirs(os.path.join(dir, "wurcs"))
except OSError:
    pass

try:
    os.makedirs(os.path.join(dir, "redend"))
except OSError:
    pass


def writter(acc, seq_type, seq):
    fn = os.path.join(dir, seq_type, acc+'.txt')
    wh = open(fn, "w")
    wh.write(seq.strip()+"\n")
    wh.close()
    return acc

seen=set()
for m in w.itermotif():
    if m.get("collection") != AllMotif.id:
        continue

    glytoucan = m.get('glytoucan')
    if glytoucan in seen:
        continue
    seen.add(glytoucan)
    wseq = m.get("wurcs")
    gseq = m.get("glycoct")
    red = m.get("redend")
    red = str(red)

    if wseq:
        writter(glytoucan, "wurcs", wseq)

    if gseq:
        writter(glytoucan, "glycoct", gseq)

    writter(glytoucan, "redend", red)

    print glytoucan
