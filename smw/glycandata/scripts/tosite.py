#!/bin/env python27

from getwiki import GlycanDataWiki, GlycanDataDiskCache

import sys

w = GlycanDataWiki()

d = GlycanDataDiskCache(sys.argv[1])
sys.argv.pop(1)

fr = None; to = None
while len(sys.argv) > 1:
    if sys.argv[1] == "--from":
	fr = sys.argv[2]
	sys.argv.pop(1); sys.argv.pop(1)
    elif sys.argv[1] == "--to":
	to = sys.argv[2]
	sys.argv.pop(1); sys.argv.pop(1)

d.towiki(w,fr=fr,to=to)
