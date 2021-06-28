#!/bin/env python2

from getwiki import GPTWiki
import sys
w = GPTWiki()
w.dumpsite(sys.argv[1],exclude_categories=['Peptide','ProteinSite','Transition','TransitionGroup'])
