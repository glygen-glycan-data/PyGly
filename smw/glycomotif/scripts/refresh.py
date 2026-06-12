#!/bin/env python3.12

import sys
from getwiki import GlycoMotifWiki, GlyTouCanMotif

w = GlycoMotifWiki()

if len(sys.argv) > 1:

  if sys.argv[1] == "-":

    for p in w.iterpages(exclude_categories=['Motif','Collection']):
      print(p.name,file=sys.stderr)
      w.refresh(p)

  elif sys.argv[1] == "Publications":

    for p in w.iterpages(include_categories=['Publication']):
      print(p.name,file=sys.stderr)
      w.refresh(p)

  elif sys.argv[1] == "Enzymes":

    for p in w.iterpages(include_categories=['Enzyme']):
      print(p.name,file=sys.stderr)
      w.refresh(p)

  elif sys.argv[1] == "Collections":

    for p in w.iterpages(include_categories=['Collection']):
      print(p.name,file=sys.stderr)
      w.refresh(p)

  else:

    for p in w.iterpages(regex=sys.argv[1]):
      print(p.name,file=sys.stderr)
      w.refresh(p)

else:

  for p in w.iterpages(include_categories=['Motif']):
    print(p.name,file=sys.stderr)
    w.refresh(p)

  for p in w.iterpages(include_categories=['Collection']):
    print(p.name,file=sys.stderr)
    w.refresh(p)

