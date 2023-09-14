#!/bin/env python3

import sys, re
import findpygly
from pygly.GlycanResource import GlyTouCan

gtc = GlyTouCan(usecache=True,prefetch=False)
for acc in sys.argv[1:]:
    seq = gtc.getseq(acc,'wurcs')
    seq = re.sub(r'-([12])[ab]',r'-\1x',seq)
    seq = re.sub(r'([a-z])\d-([a-z]\d)',r'\1?-\2',seq)
    print(seq)
