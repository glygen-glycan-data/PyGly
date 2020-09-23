#!/bin/env python27

from getwiki import GPTWiki, Protein

import sys, urllib, string, csv, os.path, json
from collections import defaultdict
import Bio.SeqIO
from util import peptide_mw, mod_mw

w = GPTWiki()

spectra = sys.argv[1]
sample = sys.argv[2]
method = sys.argv[3]
w.addacquisition(name=spectra,method=method,sample=sample,type="DIA")

for jsonfile in sys.argv[4:]:
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

  tg,mod = w.addtransgroup(peptide=pepid,z1=z1,spectra=spectra,mz1=mz1,nrt=nrt,rt=rt,prt=rt,
                           transitions=trans,ntransition=len(trans))
