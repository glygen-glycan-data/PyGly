#!/bin/env python3.12

from getwiki import GlycanDataWiki, GlycanDataDiskCache

import sys

w = GlycanDataWiki()

d = GlycanDataDiskCache(sys.argv[1])
d.tocache(w)
