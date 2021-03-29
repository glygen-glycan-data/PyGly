#!/bin/env python27

import sys
from getwiki import GPTWiki

w = GPTWiki()

if len(sys.argv) > 1:

  if sys.argv[1] == "-":

    for p in w.iterpages(exclude_categories=['Transition','Peptide','Protein','TransitionGroup','Glycan']):
      print >>sys.stderr, p.name
      w.refresh(p)

  elif sys.argv[1] == "Peptide":
    for p in w.iterpeptideids():
      print >>sys.stderr, p
      w.refresh(p)

  elif sys.argv[1] == "Transition":
    for p in w.itertransitionids():
      print >>sys.stderr, p
      w.refresh(p)

  elif sys.argv[1] == "TransitionGroup":
    for p in w.itertransgroupids():
      print >>sys.stderr, p
      w.refresh(p)

  elif sys.argv[1] in ('Protein','Glycan'):

    for p in w.iterpages(include_categories=sys.argv[1:]):
      print >>sys.stderr, p.name
      w.refresh(p)

  else:

    for p in w.iterpages(regex=sys.argv[1]):
      print >>sys.stderr, p.name
      w.refresh(p)

else:

  for p in w.iterpages(include_categories=('Transition','Peptide','Protein','TransitionGroup','Glycan')):
      print >>sys.stderr, p.name
      w.refresh(p)

