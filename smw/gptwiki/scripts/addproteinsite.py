#!/bin/env python2

from getwiki import GPTWiki, Peptide, ProteinSite

import sys

w = GPTWiki()
for pep in w.iterpep():
    for al in pep.get('alignments',[]):
	pr = al.get('protein')
	prsites = al.get('prsites',"").split('|')
	for prs in prsites:
	    aa = prs[0]
            pos = int(prs[1:])
            ps = ProteinSite(protein=pr,aa=aa,position=pos)
            if w.put(ps):
		print ProteinSite.pagename(protein=pr,aa=aa,position=pos)
	    al.append('site',ps)
	# al.delete('prsites')
    if w.put(pep):
	print pep.get('id')
