#!/bin/env python3.12

import sys, csv
from collections import defaultdict 

from getwiki import GlycanData, Glycan
w = GlycanData()

def read_gtc_taxid(fn):
    taxa = defaultdict(list)
    reader = csv.DictReader(open(fn),dialect='excel-tab')
    for r in reader:
        if r['Source'] != "GlyTouCan":
            continue
        acc = r['GlyTouCanAccession']
        tid = r['TaxID']
        taxa[acc].append(tid)
    return taxa

taxa = read_gtc_taxid(sys.argv[1])
    
for m in w.iterglycan():
    acc = m.get('accession')
    for ta in taxa[acc]:
        m.add_annotation(value=ta,property="Taxonomy",source="GlyTouCan",type="Taxonomy",sourceid=acc)
    if w.put(m):
        print(acc)
