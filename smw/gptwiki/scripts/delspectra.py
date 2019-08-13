#!/bin/env python27

from getwiki import GPTWiki, Peptide
import sys
from collections import defaultdict

w = GPTWiki()
for tgpage in w.iterpages(include_categories=['TransitionGroup']):
    tg = w.get(tgpage.name)
    if tg.get('spectra') in sys.argv[1:]:
	print >>sys.stderr, "Delete transition group",tgpage.name
	w.delete(tgpage.name)
for sp in sys.argv[1:]:
    print >>sys.stderr, "Delete spectra",sp
    w.delete(sp)
