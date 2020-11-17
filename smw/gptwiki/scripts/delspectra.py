#!/bin/env python27

from getwiki import GPTWiki, Peptide
import sys
from collections import defaultdict

w = GPTWiki()
for sp in sys.argv[1:]:
  for tg in w.iterspectgs(sp):
    if not tg:
      continue
    print >>sys.stderr, "Delete transition group",tg.get('id')
    w.delete(tg.get('id'))
  print >>sys.stderr, "Delete spectra",sp
  w.delete(sp)
