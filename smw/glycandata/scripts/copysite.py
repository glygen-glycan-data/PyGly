#!/bin/env python3.12

from getwiki import GlycanDataWiki
import sys

fromwiki = sys.argv[1].upper()
towiki = sys.argv[2].upper()

assert fromwiki in ("PROD","DEV","TEST")
assert towiki in ("PROD","DEV","TEST")
assert fromwiki != towiki

w1 = GlycanDataWiki(smwenv=fromwiki, quiet=True)
print >>sys.stderr, "from: %s"%(w1.title(),)
w2 = GlycanDataWiki(smwenv=towiki, quiet=True)
print >>sys.stderr, "  to: %s"%(w2.title(),)

dummy = raw_input("Enter to proceed, <Ctrl-C> to abort:")

currentglycan = set()
for m1 in w1.iterglycan():
    id = m1.get('id')
    assert(id)
    currentglycan.add(id)
    if w2.put(m1):
        print >>sys.stderr, "Pushing %s to %s"%(id,w2.title())
    else:
        print >>sys.stderr, "No change to %s in %s"%(id,w2.title())

for m2 in w2.iterglycan():
    id = m2.get('id')
    assert(id)
    if id not in currentglycan:
        w2.delete(id)
        print >>sys.stderr, "Delete %s from %s"%(id,w2.title())

for an in w2.iterannotation():
    id = an.get('id')
    if an.get('hasglycan') not in currentglycan:
        w2.delete(id)
        print >>sys.stderr, "Delete %s from %s"%(id,w2.title())
