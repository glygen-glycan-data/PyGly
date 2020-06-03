#!/bin/env python27

import sys, time, traceback
from collections import defaultdict

from getwiki import GlycanData
w = GlycanData()

import findpygly
from pygly.GlycanResource import GlycoMotif

gm = GlycoMotif(prefetch=True,verbose=False)

collections = ('GGM','GTC')

allmotifs = dict()
for coll in collections:
    for acc,gtcacc,redend,aglycon,names in gm.allmotifs(coll):
	if aglycon == None:
	    aglycon = ""
        allmotifs[acc] = dict(names=names,aglycon=aglycon,redend=redend,gtcacc=gtcacc)

for g in w.iterglycan():

    start = time.time()
    acc = g.get('accession')

    g.delete_annotations(type='Motif')
    motifs = set()
    for coll in collections:
        motifs.update(gm.getmotif(coll,acc))
    g.set_annotation(value=sorted(motifs),
                     property='Motif',
                     source='GlycoMotif',
                     type='Motif')
    namedmotifs = set()
    for m in motifs:
        namedmotif = ":".join([m, 
                               '//'.join(allmotifs[m]['names']),
                               allmotifs[m]['aglycon'],
                               allmotifs[m]['redend'],
                               allmotifs[m]['gtcacc']])
	namedmotifs.add(namedmotif)
        
    g.set_annotation(value = sorted(namedmotifs),
                     property='NamedMotif',
                     source='GlycoMotif',
                     type='Motif')
                                    
    if w.put(g):
        print >>sys.stderr, "%s updated in %.2f sec"%(g.get('accession'),time.time()-start,)
    else:
	print >>sys.stderr, "%s checked in %.2f sec"%(g.get('accession'),time.time()-start,)

