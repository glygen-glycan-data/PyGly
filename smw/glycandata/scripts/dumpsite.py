#!/bin/env python3.12

from getwiki import GlycanData
import sys
w = GlycanData()
if sys.argv[1] == "--all":
    sys.argv.pop(1)
    w.dumpsite(sys.argv[1])
else:
    w.dumpsite(sys.argv[1],exclude_categories=['Glycan','Annotation'])
