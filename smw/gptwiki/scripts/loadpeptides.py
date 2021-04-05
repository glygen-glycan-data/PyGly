#!/bin/env python27

from getwiki import GPTWiki, Alignment

import sys, urllib, string, csv
import Bio.SeqIO
from util import peptide_mw, mod_mw
from operator import itemgetter

w = GPTWiki()
gmw = dict()
gsym = dict()
seen = set()
for peptidefile in sys.argv[1:]:
  rest,sample,method,index,extn = peptidefile.rsplit('.',4)
  for l in csv.DictReader(open(peptidefile),dialect='excel-tab'):

    seq = l['PeptideSequence']
    glyspec = l['Glycans']
    modspec = l['Mods']
    praccs = l['ProteinName']
    pepid = l.get('PeptideID')

    if (seq,glyspec,modspec) in seen:
	continue

    if '?' in glyspec:
	continue

    seen.add((seq,glyspec,modspec))
    # print >>sys.stderr, seq,glyspec,modspec

    glycan = []
    if glyspec != "-":
        glycan = map(lambda t: (t[1],seq[int(t[0])-1]+str(t[0])),map(lambda s: s.split(':'),glyspec.split(',')))

    badglys = set()
    for glyacc in set(map(itemgetter(0),glycan)):
	if not w.get(glyacc):
	    badglys.add(glyacc)

    if len(badglys)>0:
	print >>sys.stderr, "Warning: Can't resolve glycan accession(s):",", ".join(sorted(badglys))
	continue

    alignments = []
    for pracc in map(str.strip,praccs.split(',')):

      if not w.get(pracc):
	print >>sys.stderr, "Warning: Can't resolve protein accession:",pracc
        continue

      prseq = "".join(w.get(pracc).get('sequence').split())
      pos = 0
      while pos < len(seq):
        try:
            st = prseq.index(seq,pos)
	    prsites = "|".join(map(lambda t: t[1][0]+str(st+int(t[1][1:])),glycan))
            alignments.append(Alignment(protein=pracc,start=st+1,end=st+len(seq),laa=prseq[st-1] if st > 0 else "-",
                                        raa=prseq[st+len(seq)] if (st+len(seq)) < len(prseq) else "-",prsites=prsites))
	    pos = st+1
        except ValueError:
	    break
    if len(alignments) == 0:
        print >>sys.stderr, "Warning: Can't find peptide %s in protein's sequence for accession: %s"%(seq,pracc)
        continue

    laa = alignments[0].get('laa')
    raa = alignments[0].get('raa')
    
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
    if modspec not in ("-",""):
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
    name = "%s.%s.%s"%(laa,name,raa)
    # for i,m in enumerate(glycan):
    #	name += "%s%s@%s"%(("+" if i == 0 else ","),m[0],m[1])
    nox = name.count('[Ox]')
    if pepid:
      p,modified = w.addpeptide(sequence=seq,glycans=glycan,mods=mods,mw=mw,name=name,alignments=alignments,nox=nox,id=pepid)
    else:
      p,modified = w.addpeptide(sequence=seq,glycans=glycan,mods=mods,mw=mw,name=name,alignments=alignments,nox=nox)
    if modified:
        print p.get('id')
