#!/bin/env python27

import findpygly
from pygly.GlycoMotifWiki import GlycoMotifWiki

w = GlycoMotifWiki()
w.deletemany(category='Motif',verbose=True)
