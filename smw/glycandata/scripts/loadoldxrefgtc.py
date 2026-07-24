#!/bin/env python3.12

import sys, csv
from collections import defaultdict 

from getwiki import GlycanData, Glycan
w = GlycanData()

def read_crossref(fn,prop):
    headers = None
    xrefs = defaultdict(set)
    glyconnect = False
    glycomedb = False
    reader = csv.DictReader(open(fn),dialect='excel-tab')
    for r in reader:
        if not headers:
            headers = reader.fieldnames
            if len(headers) > 2:
                if headers[2] == "GlyConnectAccessionType":
                    glyconnect = True
                elif headers[2] == "GlycOID":
                    glycomedb = True
        # rule for GlyConnect
        if glyconnect and r[headers[2]] != prop:
            continue
        # rule for GlycomeDB and GlycO
        if glycomedb and prop == "GlycO":
            xref = r[headers[2]]
        else:
            xref = r[headers[1]]
        if not xref.strip():
            continue
        acc = r['GlyTouCanAccession']
        xrefs[acc].add(xref)
    return xrefs

prop = sys.argv[1]
xrefs = read_crossref(sys.argv[2],prop)
    
for m in w.iterglycan():
    acc = m.get('accession')
    m.delete_annotations(property=prop,type="CrossReference",source="GlyTouCan")
    for xr in xrefs[acc]:
        m.add_annotation(value=xr,property=prop,source="EdwardsLab",type="CrossReference")
    if w.put(m):
        print(acc)
