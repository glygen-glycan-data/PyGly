#!/bin/env python3.12

import sys, csv
from collections import defaultdict 

from getwiki import GlycanData, Glycan
w = GlycanData()

def read_crossref(fn):
    headers = None
    xrefs = defaultdict(set)
    reader = csv.DictReader(open(fn),dialect='excel-tab')
    for r in reader:
        if not headers:
            headers = reader.fieldnames
        acc = r['GlyTouCanAccession']
	xref = r[headers[1]]
        xrefs[acc].add(xref)
    return xrefs

prop = sys.argv[1]
xrefs = read_crossref(sys.argv[2])
    
for m in w.iterglycan():
    acc = m.get('accession')
    m.delete_annotations(property=prop,type="CrossReference",source="GlyTouCan")
    for xr in xrefs[acc]:
        m.add_annotation(value=xr,property=prop,source="GlyTouCan",type="CrossReference")
    if w.put(m):
        print(acc)
