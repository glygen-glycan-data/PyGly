#!/bin/env python27

import sys

from getwiki import GlycanData
w = GlycanData()

def accessions(args):
    if len(args) == 0:
        for it in sys.stdin:
            yield it.strip()
    else:
        for fn in args:
            for it in open(fn):
                yield it.strip()

current_glygen = set(accessions(sys.argv[1:]))

for acc in w.iterglycanid():
    m = w.get(acc)
    if acc in current_glygen:
        m.set_annotation(value=acc,property="GlyGen",source="EdwardsLab",type="CrossReference")
    else:
        m.delete_annotations(property="GlyGen",source="EdwardsLab",type="CrossReference")
    if w.put(m):
        print acc,"updated"
    else:
        print acc,"checked"
