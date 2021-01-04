#!/bin/env python27

from getwiki import GlycoMotifWiki
import sys

w = GlycoMotifWiki()

m = w.get(sys.argv[1])
if m:
    coll,acc = sys.argv[2].split('.')
    m.delete('sameas')
    m.set('accession',acc)
    m.set('collection',coll)
    m.set('id',sys.argv[2])
    if w.put(m):
	print "Copied:",sys.argv[1],"to",m.get('id')
