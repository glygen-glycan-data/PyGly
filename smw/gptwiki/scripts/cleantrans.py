#!/bin/env python27

from getwiki import GPTWiki, Peptide
import sys
from collections import defaultdict

w = GPTWiki()
tg2trans = defaultdict(set)
peptranscnt = defaultdict(int)
currtrans = set()
for tgpage in w.iterpages(include_categories=['TransitionGroup']):
    tg = w.get(tgpage.name)
    tgid = tg.get('id')
    for t,i in tg.get('transitions'):
	tr = w.get(t)
	if tr:
	    tg2trans[tgid].add(t)
	    currtrans.add(t)
    pep = tg.get('peptide')
    peptranscnt[pep] += 1
for tgid in tg2trans:
    if len(tg2trans[tgid]) < 4:
	print >>sys.stderr, "Delete transition group",tgid
	w.delete(tgid)
for peppage in w.iterpages(include_categories=['Peptide']):
    pep = w.get(peppage.name)
    pepid = pep.get('id')
    if peptranscnt[pepid] == 0:
	print >>sys.stderr, "Delete peptide "+pepid
        w.delete(pepid)
for tpage in w.iterpages(include_categories=['Transition']):
    t = w.get(tpage.name)
    tid = t.get('id')
    if tid not in currtrans:
	print >>sys.stderr, "Delete transition "+tid
	w.delete(tid)
