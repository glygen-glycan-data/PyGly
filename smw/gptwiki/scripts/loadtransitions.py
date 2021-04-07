#!/bin/env python27

from getwiki import GPTWiki, Protein

import sys, urllib, string, csv, os.path
from collections import defaultdict
import Bio.SeqIO
from util import peptide_mw, mod_mw

def asscan(s):
    t = s.rstrip(')').split('(')
    return [int(t[0])]+t[1].split(',')

w = GPTWiki()

spectra2tg = defaultdict(set)
for tgpage in w.iterpages(include_categories=['TransitionGroup']):
    tg = w.get(tgpage.name)
    spectra = tg.get('spectra')
    spectra2tg[spectra].add(tg.get('id'))

allspectra = set()

for transfile in sys.argv[1:]:
  spectra,sample,method,index,extn = transfile.rsplit('.',4)
  spectra = os.path.split(spectra)[1]
  allspectra.add(spectra)
  w.addacquisition(name=spectra,method=method,sample=sample)
  tgroup = defaultdict(dict)
  for l in csv.DictReader(open(transfile),dialect='excel-tab'):
    gphash = l.get('GlycopeptideHash')
    seq = l['PeptideSequence']
    glyspec = l['Glycans']
    modspec = l['Mods']
    glycan = []
    if '?' in glyspec:
	continue
    if glyspec != "-":
        glycan = map(lambda t: (t[1],seq[int(t[0])-1]+str(t[0])),map(lambda s: s.split(':'),glyspec.split(',')))
    mods = []
    if modspec != "-":
	for mstr in modspec.split(','):
	    pos, delta = mstr.split(':')
	    delta, site, sym = mod_mw(float(delta),int(pos),seq)
	    mods.append((delta,site))
    key,pid = w.findpeptide(sequence=seq,glycans=glycan,mods=mods)
    if not pid:
	print >>sys.stderr, "Warning: Can't resolve peptide %s, %s, %s"%(seq,glycan,mods)
	continue
    rt = float(l['RetentionTime'])
    nrt = l['NormalizedRetentionTime']
    if nrt:
	nrt = float(nrt)
    prt = l.get('PeakRT')
    if prt:
	prt = float(prt)
    mz1 = float(l['PrecursorMz'])
    z1 = int(l['PrecursorCharge'])
    mz2 = float(l['ProductMz'])
    z2 = int(l['ProductCharge'])
    label = l['Annotation']
    relint = float(l['LibraryIntensity'])
    label = label.rstrip('+')
    scans = map(asscan,l['Scans'].split(';'))

    t,modified = w.addtransition(peptide=pid,label=label,mz1=mz1,z1=z1,mz2=mz2,z2=z2)
    if modified:
        print t.get('id')

    if (pid,z1,spectra) not in tgroup:
        tgroup[(pid,z1,spectra)] = dict(transitions=[],nrt=nrt,rt=rt,prt=prt,mz1=mz1,scans=scans,gphash=gphash)
    tgroup[(pid,z1,spectra)]['transitions'].append((t.get('id'),relint))
    
  for pid,z1,spectra in tgroup:
    tgroup[(pid,z1,spectra)]['ntransition'] = len(tgroup[(pid,z1,spectra)]['transitions'])
    tg,mod = w.addtransgroup(peptide=pid,z1=z1,spectra=spectra,**tgroup[(pid,z1,spectra)])
    if mod:
        print tg.get('id')
    if tg.get('id') in spectra2tg[spectra]:
        spectra2tg[spectra].remove(tg.get('id'))

for spectra in allspectra:
  for tgid in spectra2tg[spectra]:
    if w.cleartransgroup(tg):                                                                                              
      print "Clear",tg.get("id")                                                                                           
