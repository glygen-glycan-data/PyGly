#!/bin/env python27

import sys
from getwiki import GPTWiki

w = GPTWiki()

if len(sys.argv) > 1:

  if sys.argv[1] == "-":

    for p in w.iterpages(exclude_categories=['Transition','Peptide','Protein','TransitionGroup','Glycan']):
      print >>sys.stderr, p.name
      w.refresh(p)

  elif sys.argv[1] in ('Transition','Peptide','Protein','Glycan','TransitionGroup'):

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

