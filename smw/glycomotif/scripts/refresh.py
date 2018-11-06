#!/bin/env python27

import sys
from getwiki import GlycoMotifWiki, GlyTouCanMotif

w = GlycoMotifWiki()

if len(sys.argv) > 1:

  for p in w.iterpages(regex=sys.argv[1]):
    print >>sys.stderr, p.name
    w.refresh(p)

else:

  for p in w.iterpages(include_categories=['Motif']):
    print >>sys.stderr, p.name
    w.refresh(p)

  for p in w.iterpages(include_categories=['Collection']):
    print >>sys.stderr, p.name
    w.refresh(p)

