#!/bin/env python2

import sys, time, traceback
from collections import defaultdict

from getwiki import GlycanData
w = GlycanData()

import findpygly
from pygly.GlycanResource import GlycoMotif

gm = GlycoMotif(prefetch=True,verbose=False,usecache=False)

collections = ('GGM',)

allmotifs = dict()
for coll in collections:
    for acc,gtcacc,alignment,redend,aglycon,names,pmids,keywords,dbxrefs in gm.allmotifs(coll):
	if aglycon == None:
	    aglycon = ""
	if coll == "GGM" and int(acc.split('.',1)[1]) > 1000:
	    classification = True
	else:
	    classification = False
        allmotifs[acc] = dict(names=names,alignment=alignment,aglycon=aglycon,redend=redend,gtcacc=gtcacc,pmids=pmids,keywords=keywords,classification=classification,dbxrefs=dbxrefs)

for g in w.iterglycan():

    start = time.time()
    acc = g.get('accession')

    g.delete_annotations(type='Motif')
    strictglygenmotifs = set()
    allglygenmotifs = set()
    classmotifs = set()
    strictclassmotifs = set()
    names = set()
    for coll in collections:
	for m,s in gm.getmotif(coll,acc):
	    if not allmotifs[m]['classification']:
		if s:
		    strictglygenmotifs.add(m)
		    if allmotifs[m]['alignment'] == "Whole-Glycan":
			names.add((m,allmotifs[m]['names'][0]))
	        allglygenmotifs.add(m)
	    else:
		if s:
		    strictclassmotifs.add(m)
                classmotifs.add(m)
    g.set_annotation(value=sorted(strictglygenmotifs),
                     property='StrictMotif',
                     source='GlycoMotif',
                     type='Motif')
    g.set_annotation(value=sorted(allglygenmotifs),
                     property='Motif',
                     source='GlycoMotif',
                     type='Motif')
    g.set_annotation(value=sorted(classmotifs),
                     property='ClassMotif',
                     source='GlycoMotif',
                     type='Motif')
    g.set_annotation(value=sorted(strictclassmotifs),
                     property='StrictClassMotif',
                     source='GlycoMotif',
                     type='Motif')

    for m,n in names:
        g.set_annotation(value=n,property="SemanticName",source="GlycoMotif",sourceid=m,type="Name")

    # namedmotifs = set()
    # for m in motifs:
    #     namedmotif = ":".join([m, 
    #                            '//'.join(allmotifs[m]['names']),
    #                            allmotifs[m]['aglycon'],
    #                            allmotifs[m]['redend'],
    #                            allmotifs[m]['gtcacc']])
    #	namedmotifs.add(namedmotif)
    #    
    # g.set_annotation(value = sorted(namedmotifs),
    #                  property='NamedMotif',
    #                  source='GlycoMotif',
    #                  type='Motif')
                                    
    if w.put(g):
        print >>sys.stderr, "%s updated in %.2f sec"%(g.get('accession'),time.time()-start,)
    else:
	print >>sys.stderr, "%s checked in %.2f sec"%(g.get('accession'),time.time()-start,)

