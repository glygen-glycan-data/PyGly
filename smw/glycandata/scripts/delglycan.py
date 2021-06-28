#!/bin/env python2

import re, sys
from getwiki import GlycanData

w = GlycanData()
for l in sys.stdin:
    print l.strip()
    try:
        w.delete(l.strip())
    except OSError:
	pass
