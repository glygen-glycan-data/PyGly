#!/bin/env python27

import sys

from getwiki import GlycanDataWiki, Glycan
w = GlycanDataWiki()

def accessions(args):
    if len(args) == 0:
        for it in sys.stdin:
            yield it.strip()
    else:
        for fn in args:
            for it in open(fn):
                yield it.strip()

current_glygen = set(accessions(sys.argv[1:]))

for m in w.iterglycan():
    m.delete_annotations(property="GlyGen",source="EdwardsLab",type="CrossReference")
    acc = m.get('accession')
    if acc in current_glygen:
        m.set_annotation(value=acc,property="GlyGen",source="EdwardsLab",type="CrossReference")
    if w.put(m):
        print acc
