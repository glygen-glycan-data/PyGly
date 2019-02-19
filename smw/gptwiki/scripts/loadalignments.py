#!/bin/env python27

from getwiki import GPTWiki, Protein

import sys, urllib, string
from collections import defaultdict
import Bio.SeqIO

w = GPTWiki()
alignfile = sys.argv[1]
alignments = defaultdict(list)
for l in open(alignfile):
    sl = l.split()
    st = int(sl[1])+1
    ed = int(sl[2])
    pep = sl[4]
    laa = sl[3]
    raa = sl[5]
    pracc = sl[12][1:]
    alignments[pep].append((pracc,st,ed))

for p in w.iterpeptides():
    seq = p.get('sequence')
    if seq in alignments:
	p.update(alignment=alignments[seq])
	if w.put(p):
	    print >>sys.stderr, p.get('id')
