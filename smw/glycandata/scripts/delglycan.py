#!/bin/env python27

import re, sys
from getwiki import GlycanData

w = GlycanData()
for l in sys.stdin:
    print l.strip()
    w.delete(l.strip())
