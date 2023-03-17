#!/bin/env python

import sys, os, os.path
import findpygly
from pygly.GlycanFormatter import WURCS20Format, GlycanParseError

fmt = WURCS20Format()

pos = int(sys.argv[1])

for l in sys.stdin:
    wurcs = l.split()[pos]
    try:
        g = fmt.toGlycan(wurcs)
        if not g.repeated():
            comp = g.native_elemental_composition()
            comp['H'] += 2
            comp['O'] += 1
            print("%s\t%s"%(l.strip(),comp.compactstr()))
    except (GlycanParseError,KeyError):
        pass
