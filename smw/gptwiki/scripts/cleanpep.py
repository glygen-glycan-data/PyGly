#!/bin/env python27

from getwiki import GPTWiki, Peptide

import sys, urllib, string
import Bio.SeqIO
from util import peptide_mw, mod_mw

seen = set()

w = GPTWiki()
peps = []
for p in w.iterpages(include_categories=['Peptide']):
    pep = w.get(p.name)
    pepkey = Peptide.key(pep.get('sequence'),pep.get('glycan',[]),pep.get('mod',[]))
    if pepkey in seen:
	print >>sys.stderr, p.name
	w.delete(p.name)
    seen.add(pepkey)
