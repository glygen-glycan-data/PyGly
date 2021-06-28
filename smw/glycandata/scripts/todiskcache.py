#!/bin/env python2

from getwiki import GlycanDataWiki, GlycanDataDiskCache

import sys

w = GlycanDataWiki()

d = GlycanDataDiskCache(sys.argv[1])
d.tocache(w)
