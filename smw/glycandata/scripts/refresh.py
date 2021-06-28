#!/bin/env python2

import sys
from getwiki import GlycanData

w = GlycanData()

if len(sys.argv) > 1:

  if sys.argv[1] == "-":

    for p in w.iterpages(exclude_categories=['Glycan','Annotation']):
      print >>sys.stderr, p.name
      w.refresh(p)

  elif sys.argv[1] == "stdin":

    for p in map(str.strip,sys.stdin):
      print >>sys.stderr, p
      w.refresh(p)

  else:

    for p in w.iterpages(regex=sys.argv[1]):
      print >>sys.stderr, p.name
      w.refresh(p)

else:

  for p in w.iterpages(include_categories=['Annotation']):
    print >>sys.stderr, p.name
    w.refresh(p)
  for p in w.iterpages(include_categories=['Glycan']):
    print >>sys.stderr, p.name
    w.refresh(p)

