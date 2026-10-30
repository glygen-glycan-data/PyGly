#!/bin/env python3.12

from getwiki import GlycoMotifWiki
import sys

w = GlycoMotifWiki()
w.loadsite(sys.argv[1])
