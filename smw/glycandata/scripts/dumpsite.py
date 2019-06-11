#!/bin/env python27

from getwiki import GlycanDataWiki
import sys
w = GlycanDataWiki()
if sys.argv[1] == "--all":
    sys.argv.pop(1)
    w.dumpsite(sys.argv[1])
else:
    w.dumpsite(sys.argv[1],exclude_categories=['Glycan','Annotation'])
