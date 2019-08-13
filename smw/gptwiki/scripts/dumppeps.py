#!/bin/env python27

from getwiki import GPTWiki, Peptide

import sys, urllib, string
import Bio.SeqIO

print "\t".join(map(str,["Spectra","Accession","Peptide","Site","Glycan","TGAccession","Charge","PrecursorMZ"]))
seen = set()
w = GPTWiki()
for tgpage in w.iterpages(include_categories=['TransitionGroup']):
    tg = w.get(tgpage.name)
    pep = w.get(tg.get('peptide'))
    pepid = pep.get('id')
    z1 = tg.get('z1')
    mz1 = tg.get('mz1')
    spectra = tg.get('spectra')
    if (spectra,pepid,z1) not in seen:
	pepseq = list(pep.get('sequence'))
	for deltastr,pos in pep.get('mod',[]):
	    aa = pos[0]
	    pos = int(pos[1:])-1
	    if round(deltastr,3) == 57.021:
		pepseq[pos] += ":m"
	    elif round(deltastr,3) in (15.995,15.996):
		pepseq[pos] += ":o"
	pepseq = "".join(pepseq)
	sites = []
	glycans = []
	for glyacc,pos in pep.get('glycan',[]):
	    aa = pos[0]
	    pos = int(pos[1:])-1
	    sites.append(pos)
	    glycans.append(glyacc)
	sites = ",".join(map(str,sites))
	glycans = ",".join(map(str,glycans))
        print "\t".join(map(str,[spectra,pepid,pepseq,sites,glycans,tg.get('id'),z1,mz1]))
	seen.add((spectra,pepid,z1))
