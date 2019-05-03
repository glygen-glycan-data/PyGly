#!/bin/env python27

from getwiki import GlycanDataWiki
import sys
w = GlycanDataWiki()
w.dumpsite(sys.argv[1],exclude_categories=['Glycan','Annotation'])
