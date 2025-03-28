#!/bin/env python3.12

import sys, time
from collections import defaultdict
import csv

from getwiki import GlycanData, Glycan
w = GlycanData()

gdb2gog = defaultdict(set)
for l in open(sys.argv[1]):
    sl = l.split()
    gdb2gog[sl[0]].add(sl[1])

for m in w.iterglycan():
    start = time.time()
    gog = set()
    gdbset = set()
    for ann in m.annotations(property="GlycomeDB"):
        for gdb in ann.get('value',[]):
            gdbset.add(gdb) 
            gog.update(gdb2gog[gdb])
    if len(gog) > 0:
        m.set_annotation(value=list(gog), property="GlycO",source="GlycomeDB",type="CrossReference")
    else:
        m.delete_annotations(source="GlycomeDB",property="GlycO",type="CrossReference")
    if w.put(m):
        print("%s updated in %.2f sec"%(m.get('accession'),time.time()-start,),file=sys.stderr)
    else:
        print("%s checked in %.2f sec"%(m.get('accession'),time.time()-start,),file=sys.stderr)
