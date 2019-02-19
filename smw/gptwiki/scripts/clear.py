#!/bin/env python27

from getwiki import GPTWiki
import sys

w = GPTWiki()
for cat in sys.argv[1:]:
    assert cat in ("Transition","TransitionGroup","Peptide","Protein","Glycan")
    w.deletemany(category=cat,verbose=True)

