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
            comp = g.iupac_composition(aggregate_basecomposition=False)
            count = comp['Count']
            del comp['Count']
            if comp['Xxx'] == 0:
                print("%s\t%s\t%s"%(l.strip(),count,comp))
    except (GlycanParseError,KeyError):
        pass
