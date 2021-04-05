#!/bin/env python2

from getwiki import GlycoMotifWiki
import sys
w = GlycoMotifWiki()
w.dumpsite(sys.argv[1],exclude_regex='^GM\.')
