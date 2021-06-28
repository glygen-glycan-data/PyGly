#!/bin/env python2

import sys, time, traceback, hashlib, os, os.path, glob
from collections import defaultdict

import findpygly
from getwiki import GlycanData
w = GlycanData()

allfn = set(glob.glob(sys.argv[1]+"/G*.txt"))
try:
    os.makedirs(sys.argv[1])
except OSError:
    pass

for g in w.iterglycan():
    acc = g.get('accession')

    seqtype = sys.argv[1].split('/')[-1]
    assert(seqtype in ('wurcs','glycoct','genglycoct'))

    filename = sys.argv[1] + "/" + acc + ".txt"

    if filename in allfn:
        allfn.remove(filename)	
	continue

    seq = None
    if seqtype == 'wurcs' and g.has_annotations(property='WURCS',type='Sequence',source='GlyTouCan'):
	seq = g.get_annotation_value(property='WURCS',type='Sequence',source='GlyTouCan')

    if seqtype == 'glycoct' and g.has_annotations(property='GlycoCT',type='Sequence',source='GlyTouCan'):
	seq = g.get_annotation_value(property='GlycoCT',type='Sequence',source='GlyTouCan')

    if seqtype == 'genglycoct' and g.has_annotations(property='GlycoCT',type='Sequence',source='EdwardsLab'):
	seq = g.get_annotation_value(property='GlycoCT',type='Sequence',source='EdwardsLab')

    if not seq:
	continue

    print >>sys.stderr, "Write:",acc
    wh = open(filename,'w')
    wh.write(seq)
    wh.close()

for fn in allfn:
    print >>sys.stderr, "Remove:",acc
    os.unlink(fn)
