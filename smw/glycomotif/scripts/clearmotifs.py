#!/bin/env python27

from getwiki import GlycoMotifWiki
import sys

w = GlycoMotifWiki()
if len(sys.argv) > 1:
    w.deletemany(regex=sys.argv[1],verbose=True)
else:
    w.deletemany(category='Motif')
