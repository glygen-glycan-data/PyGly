#!/bin/env python27

from getwiki import GlycanDataWiki
import sys

w = GlycanDataWiki()
w.loadsite(sys.argv[1])
