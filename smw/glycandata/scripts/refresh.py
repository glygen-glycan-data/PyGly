#!/bin/env python3.12

import sys
from getwiki import GlycanData

w = GlycanData()

if len(sys.argv) > 1:

  if sys.argv[1] == "-":

    for p in w.iterpages(exclude_categories=['Glycan']):
      print(p.name,file=sys.stderr)
      w.refresh(p)

  elif sys.argv[1] == "stdin":

    for p in map(str.strip,sys.stdin):
      print(p,file=sys.stderr)
      w.refresh(p)

  else:

    for p in w.iterpages(regex=sys.argv[1]):
      print(p.name,file=sys.stderr)
      w.refresh(p)

else:

  for p in w.iterpages(include_categories=['Glycan']):
    print(p.name,file=sys.stderr)
    w.refresh(p)

