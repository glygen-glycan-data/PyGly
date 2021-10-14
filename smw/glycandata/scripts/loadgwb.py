#!/bin/env python2

import sys, time, traceback, hashlib
from collections import defaultdict

import findpygly

from pygly.GWBFormat import GWBFormat
f = GWBFormat()

from getwiki import GlycanData
w = GlycanData()

def accessions():
    if len(sys.argv) > 1:
	for arg in sys.argv[1:]:
	    g = w.get(arg.strip())
	    if g:
		yield g
    else:
	for g in w.iterglycan():
	    yield g

for g in accessions():
    start = time.time()

    glycan = g.getGlycan()

    if not glycan:
    	continue

    try:
        wurcs = g.get_annotation_value(property="WURCS",type='Sequence',source='GlyTouCan')
    except LookupError:
	continue

    if glycan.has_root() and not glycan.repeated() and not g.has_annotations(property="GlycoWorkBench",type='Sequence',source='EdwardsLab'):
        try:
            gwb = f.toStr(wurcs)
            if gwb:
                g.set_annotation(property="GlycoWorkBench",value=gwb,type='Sequence',source='EdwardsLab')
        except:
            traceback.print_exc()

    if w.put(g):
        print >>sys.stderr, "%s updated in %.2f sec"%(g.get('accession'),time.time()-start,)
    else:
	print >>sys.stderr, "%s checked in %.2f sec"%(g.get('accession'),time.time()-start,)

