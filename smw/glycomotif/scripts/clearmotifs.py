#!/bin/env python27

from getwiki import GlycoMotifWiki
import sys

w = GlycoMotifWiki()
w.deletemany(category='Motif',verbose=True)

