#!/bin/env python27

import sys, time
from collections import defaultdict
import csv

from getwiki import GlycanDataWiki, Glycan
w = GlycanDataWiki()

gdb2gog = defaultdict(set)
for l in open(sys.argv[1]):
    sl = l.split()
    gdb2gog[sl[0]].add(sl[1])

for m in w.iterglycan():
    start = time.time()
    gog = set()
    try:
        for gdb in m.get_annotation(property="GlycomeDB").get('value',[]):
	    gog.update(gdb2gog[gdb])
    except LookupError:
	pass
    if len(gog) > 0:
        m.set_annotation(value=list(gog), property="GlycO",source="GlycomeDB",type="CrossReference")
    else:
        m.delete_annotations(source="GlycomeDB",property="GlycO",type="CrossReference")
    if w.put(m):
        print >>sys.stderr, "%s updated in %.2f sec"%(m.get('accession'),time.time()-start,)
    else:
        print >>sys.stderr, "%s checked in %.2f sec"%(m.get('accession'),time.time()-start,)
