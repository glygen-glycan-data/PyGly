#!/bin/env python3.12

import re, sys
from getwiki import GlycanData

w = GlycanData()
for l in sys.stdin:
    print(l.strip())
    try:
        w.delete(l.strip())
    except OSError:
	pass
