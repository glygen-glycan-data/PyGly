#!/bin/env python27

from getwiki import GlycoMotifWiki
import sys

fromwiki = sys.argv[1].upper()
towiki = sys.argv[2].upper()

assert fromwiki in ("PROD","DEV","TEST")
assert towiki in ("PROD","DEV","TEST")
assert fromwiki != towiki

w1 = GlycoMotifWiki(smwenv=fromwiki, quiet=True)
# print >>sys.stderr, "from: %s"%(w1.title(),)
w2 = GlycoMotifWiki(smwenv=towiki, quiet=True)
# print >>sys.stderr, "  to: %s"%(w2.title(),)

import difflib

w1only = set()
w2only = set()
both = set()
for p1 in w1.iterall():
    p2 = w2._get(p1.name)
    if not p2 or not p2.exists:
	w1only.add(p1.name)
    both.add(p1.name)
    t1 = p1.text().encode('utf-8').splitlines()
    t2 = p2.text().encode('utf-8').splitlines()
    fromfile = "%s/%s"%(fromwiki,p1.name)
    tofile = "%s/%s"%(towiki,p2.name)
    for l in difflib.unified_diff(t1,t2,fromfile=fromfile,tofile=tofile):
	sys.stdout.write((l.rstrip()+'\n'))

for p2 in w2.iterall():
    if p2.name not in both:
	w2only.add(p2.name)

for n in w1only:
    print "Only in %s: %s"%(fromwiki,n)
for n in w2only:
    print "Only in %s: %s"%(towiki,n)

