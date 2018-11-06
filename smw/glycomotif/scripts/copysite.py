#!/bin/env python27

from getwiki import GlycoMotifWiki
import sys

fromwiki = sys.argv[1].upper()
towiki = sys.argv[2].upper()

assert fromwiki in ("PROD","DEV","TEST")
assert towiki in ("PROD","DEV","TEST")
assert fromwiki != towiki

w1 = GlycoMotifWiki(smwenv=fromwiki, quiet=True)
print >>sys.stderr, "from: %s"%(w1.title(),)
w2 = GlycoMotifWiki(smwenv=towiki, quiet=True)
print >>sys.stderr, "  to: %s"%(w2.title(),)

dummy = raw_input("Enter to proceed, <Ctrl-C> to abort:")

currentcoll = set()
for c1 in w1.itercollection():
    id = c1.get('id')
    assert(id)
    currentcoll.add(id)
    if w2.put(c1):
        print >>sys.stderr, "Pushing %s to %s"%(id,w2.title())
    else:
        print >>sys.stderr, "No change to %s in %s"%(id,w2.title())

currentmotif = set()
for m1 in w1.itermotif():
    id = m1.get('id')
    assert(id)
    currentmotif.add(id)
    if w2.put(m1):
        print >>sys.stderr, "Pushing %s to %s"%(id,w2.title())
    else:
        print >>sys.stderr, "No change to %s in %s"%(id,w2.title())

for c2 in w2.itercollection():
    id = c2.get('id')
    assert(id)
    if id not in currentcoll:
	w2.delete(id)
	print >>sys.stderr, "Delete %s from %s"%(id,w2.title())

for m2 in w2.itermotif():
    id = m2.get('id')
    assert(id)
    if id not in currentmotif:
	w2.delete(id)
	print >>sys.stderr, "Delete %s from %s"%(id,w2.title())
