#!/bin/env python3.12

import sys, time
from collections import defaultdict
import csv

from getwiki import GlycanData, Glycan
w = GlycanData()

headers = """
GlyTouCanAccession      GlycomeDBID     GlycOID
""".split()

print("\t".join(headers))
for acc in sorted(w.iterglycanid()):
    g = w.get(acc)
    gdb = set()
    gog = set()
    for ann in g.annotations(type="CrossReference"):
        if ann.get('property') == 'GlycomeDB':
            gdb.update(ann.get('value'))
        elif ann.get('property') == 'GlycO':
            gog.update(ann.get('value'))
    if len(gdb) == 0:
        continue
    assert len(gdb) == 1 or len(gog) <= 1
    for gdbi in sorted(gdb):
        any = False
        for gogi in sorted(gog):
            any = True
            print("\t".join([acc,gdbi,gogi]))
        if not any:
            print("\t".join([acc,gdbi,""]))
