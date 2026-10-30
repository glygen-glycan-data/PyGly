#!/bin/env python3.12

import sys, time, re
from collections import defaultdict
import csv

from getwiki import GlycanData, Glycan
w = GlycanData()

gtc2xref = defaultdict(lambda: defaultdict(set))
for f in sys.argv[1:]:
  for l in open(f):
    acc,gtcacc,source = l.split()
    gtc2xref[gtcacc][source].add(acc)

for m in w.iterglycan():
    start = time.time()
    acc = m.get('accession')

    m.delete_annotations(source="GlyGenSourceFile", type="CrossReference")
        
    for source in gtc2xref[acc]:
        if source not in ('PDB',
                          'GlyConnectStructure',
                          'GlyConnectComposition'):
            pass
        else:
            m.set_annotation(value=sorted(gtc2xref[acc][source]), property=source,
                             source="GlyGenSourceFile", type="CrossReference")

    if w.put(m):
        print("%s updated in %.2f sec"%(m.get('accession'),time.time()-start,),file=sys.stderr)
    else:
        print("%s checked in %.2f sec"%(m.get('accession'),time.time()-start,),file=sys.stderr)
