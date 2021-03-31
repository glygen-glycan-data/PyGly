#!/bin/env python27

from getwiki import GPTWiki, Protein
import getwiki

import sys, urllib, string, csv, os.path, json
from collections import defaultdict
import Bio.SeqIO
from util import peptide_mw, mod_mw

w = GPTWiki()

sample = sys.argv[1]
method = sys.argv[2]
if ":" in method:
    method,anfrac=method.split(':')
else:
    anfrac=None

tgs = defaultdict(set)
allspec = set()
lccal = defaultdict(dict)
current_spec = None
for jsonfile in sys.argv[3:]:

  spectra = os.path.split(os.path.split(jsonfile)[0])[1]
  if spectra != current_spec:
    w.addacquisition(name=spectra,method=method,anfrac=anfrac,sample=sample,xicmmu=50,type="DIA")
    allspec.add(spectra)
    current_spec = spectra

  data = json.loads(open(jsonfile).read())
  rt = data['rt']
  nrt = data.get('nrt')
  trans = []
  for s in data['series']:
      trid = s['trlib_id'].split('.')[-1]
      it = float(s['obs_relint'])/10.0
      if it <= 0:
	continue
      trans.append((trid,it))
  trans.sort(key=lambda t: -t[1])
  pepid = data['trlib_pepid']
  z1 = data['trlib_z1']
  mz1 = data['trlib_mz1']
  extras = {}
  for key in ("intensity","score","extraction","fdr"):
      if key in data:
	  extras[key] = data[key]
  tg,mod = w.addtransgroup(peptide=pepid,z1=z1,spectra=spectra,mz1=mz1,nrt=nrt,rt=rt,prt=rt,
                           transitions=trans,ntransition=len(trans),**extras)
  if 'lccalibration' in data and spectra not in lccal:
      nrtslope,nrtintercept = map(float,data['lccalibration'].split(":"))
      lccal[spectra]['nrtslope'] = nrtslope
      lccal[spectra]['nrtintercept'] = nrtintercept
  tgs[spectra].add(tg.get('id'))

for spectra in allspec:
  spec = w.get(spectra)
  if spectra in lccal:
    spec.set("nrtslope",lccal[spectra]['nrtslope'])
    spec.set("nrtintercept",lccal[spectra]['nrtintercept'])
  else:
    spec.delete("nrtslope")
    spec.delete("nrtintercept")
  w.put(spec)

  for tg in w.iterspectgs(spectra):
    if tg.get('id') not in tgs[spectra]:
      print "Removing",tg.get("id")
      # w.delete(tg.get('id'))
      tg.delete('scans')
      tg.delete('transitions')
      w.put(tg)
