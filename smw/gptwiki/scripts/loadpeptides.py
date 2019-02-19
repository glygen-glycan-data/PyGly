#!/bin/env python27

from getwiki import GPTWiki, Protein

import sys, urllib, string, csv
import Bio.SeqIO
from util import peptide_mw, mod_mw

w = GPTWiki()
gmw = dict()
gsym = dict()
seen = set()
for peptidefile in sys.argv[1:]:
  for l in csv.DictReader(open(peptidefile),dialect='excel-tab'):
    seq = l['PeptideSequence']
    glyspec = l['Glycans']
    modspec = l['Mods']
    if (seq,glyspec,modspec) in seen:
	continue
    seen.add((seq,glyspec,modspec))
    # print >>sys.stderr, seq,glyspec,modspec
    glycan = []
    if glyspec != "-":
        glycan = map(lambda t: (t[1],seq[int(t[0])-1]+str(t[0])),map(lambda s: s.split(':'),glyspec.split(',')))
    modsyms = []
    for i in range(len(glycan)):
	g = glycan[i]
	if g[0] not in gsym:
	    gi = w.get(g[0])
	    gmw[g[0]] = gi.get('mw')
	    gsym[g[0]] = gi.get('sym')
	modsyms.append(gsym[g[0]])
    mods = []
    modsmw = 0
    if modspec != "-":
	for mstr in modspec.split(','):
	    pos, delta = mstr.split(':')
	    delta, site, sym = mod_mw(float(delta),int(pos),seq)
	    mods.append((delta,site))
	    modsmw += delta
	    modsyms.append(sym)
    mw = peptide_mw(seq)+modsmw
    for g in glycan:
	mw += gmw[g[0]]
    name = list(seq)
    for m,sym in zip(glycan+mods,modsyms):
	site = int(m[1][1:])
	if sym:
	    name[site-1] = name[site-1] + "[%s]"%(sym,)
    name = "".join(name)
    # for i,m in enumerate(glycan):
    #	name += "%s%s@%s"%(("+" if i == 0 else ","),m[0],m[1])
    p,modified = w.addpeptide(sequence=seq,glycans=glycan,mods=mods,mw=mw,name=name)
    if modified:
        print p.get('id')
