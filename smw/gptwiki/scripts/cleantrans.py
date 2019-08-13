#!/bin/env python27

from getwiki import GPTWiki, Peptide
import sys
from collections import defaultdict

w = GPTWiki()
currtrans = set()
for tgpage in w.iterpages(include_categories=['TransitionGroup']):
    tg = w.get(tgpage.name)
    tgid = tg.get('id')
    for t,i in tg.get('transitions'):
	tr = w.get(t)
	if tr:
	    currtrans.add(t)
for tpage in w.iterpages(include_categories=['Transition']):
    t = w.get(tpage.name)
    tid = t.get('id')
    if tid not in currtrans:
	print >>sys.stderr, "Delete transition "+tid
	w.delete(tid)
