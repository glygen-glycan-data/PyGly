#!/bin/env python2

import sys

from getwiki import GlycanData
w = GlycanData()

for g in w.iterglycan():
    if not g.has_annotation(property="IUPAC",source="GlyTouCan",type="Sequence"):
	
    if acc in current_glygen:
        m.set_annotation(value=acc,property="GlyGen",source="EdwardsLab",type="CrossReference")
    else:
        m.delete_annotations(property="GlyGen",source="EdwardsLab",type="CrossReference")
    if w.put(m):
        print acc,"updated"
    else:
        print acc,"checked"
