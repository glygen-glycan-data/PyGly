#!/bin/env python

from getwiki import GlycoMotifWiki
import sys

w = GlycoMotifWiki()
w.loadsite(sys.argv[1])
