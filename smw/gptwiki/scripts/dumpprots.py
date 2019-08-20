#!/bin/env python27

from getwiki import GPTWiki, Peptide

import sys, urllib, string
import Bio.SeqIO

w = GPTWiki()
for pr in sorted(w.iterproteins(),key=lambda pr: pr.get('accession')):
    acc = pr.get('accession')
    desc = pr.get('description')
    sequence = "".join(pr.get('sequence').split())
    print ">%s %s"%(acc,desc)
    for i in range(0,len(sequence),60):
	print sequence[i:i+60]
