#!/bin/env python3.12

from getwiki import GlycoMotifWiki
import sys
w = GlycoMotifWiki()
w.dumpsite(sys.argv[1],exclude_regex=r'^GM\.')
