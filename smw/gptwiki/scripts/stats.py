#!/bin/env python27

from getwiki import GPTWiki, Peptide

import sys, urllib, string
import Bio.SeqIO
from collections import defaultdict

w = GPTWiki()
seenpeps = set()
sites = set()
prot2site = defaultdict(set)
samples = set()
glycans = set()
glysites = set()
site2gly = defaultdict(set)
for tg in w.itertransgroups():
    tgid = tg.get('id')
    if tg.get('peptide') in seenpeps:
	continue
    pep = w.get(tg.get('peptide'))
    if pep.get('nrt') == None:
	continue
    seenpeps.add(pep.get('id'))
    gly = pep.get('glycan')[0][0]
    glycans.add(gly)
    for al in pep.get('alignments',[]):
	site = al.get('prsites')
	prot = al.get('protein')
	print pep.get('id'),prot,site,gly
	sites.add((prot,site))
	prot2site[prot].add(site)
	site2gly[(prot,site)].add(gly)
	glysites.add((prot,site,gly))

print "Proteins:",len(prot2site)
print "Glycans:",len(glycans)
print "Sites:",len(sites)
print "GlySites:",len(glysites)
print "Glycopeptides:",len(seenpeps)
# print "Sites/Prot:",sum(len(prot2site[prot]) for prot in prot2site)/float(len(prot))
