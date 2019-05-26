#!/bin/env python27

from getwiki import GPTWiki
import sys

fromwiki = sys.argv[1].upper()
towiki = sys.argv[2].upper()

assert fromwiki in ("PROD","DEV","TEST")
assert towiki in ("PROD","DEV","TEST")
assert fromwiki != towiki

w1 = GPTWiki(smwenv=fromwiki, quiet=True)
print >>sys.stderr, "from: %s"%(w1.title(),)
w2 = GPTWiki(smwenv=towiki, quiet=True)
print >>sys.stderr, "  to: %s"%(w2.title(),)

dummy = raw_input("Enter to proceed, <Ctrl-C> to abort:")

currentids = set()
for page in w1.iterpages(include_categories=('Transition','Peptide','Protein','Glycan','TransitionGroup','Glycan')):
    id = page.name
    currentids.add(id)
    it = w1.get(id)
    if w2.put(it):
        print >>sys.stderr, "Pushing %s to %s"%(id,w2.title())
    else:
        print >>sys.stderr, "No change to %s in %s"%(id,w2.title())

for page in w2.iterpages(include_categories=('Transition','Peptide','Protein','Glycan','TransitionGroup','Glycan')):
    id = page.name
    if id not in currentids:
	w2.delete(id)
	print >>sys.stderr, "Delete %s from %s"%(id,w2.title())
