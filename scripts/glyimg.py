#!/bin/env python2
import sys, os
import findpygly
from pygly.GlycanImage import GlycanImage
from pygly.GlycanResource.GlyTouCan import GlyTouCanNoPrefetch

if len(sys.argv) <= 1:
    print >>sys.stderr, "glyimg.py [ image options ] <gtc-accession> [ <gtc-accession> ... ]"
    print >>sys.stderr, """
    Image options:
    scale        <float>                               [1.0]
    reducing_end (true|false)                          [true]
    orientation  (RL|LR|TB|BT)                         [RL]
    notation     (cfg|cfgbw|cfglink|uoxf|text|uoxfcol) [cfg]
    display      (normal|normalinfo|compact)           [normalinfo]
    format	 (png|svg)                             [png]
    """.strip()
    sys.exit(1)

gtc = GlyTouCanNoPrefetch()

imageWriter = GlycanImage()
imageWriter.verbose(True)
imageWriter.force(True)
lastopt = 0
for i in range(1,len(sys.argv),2):
    if sys.argv[i].startswith('G'):
	break
    key = sys.argv[i]
    value = sys.argv[i+1]
    imageWriter.set(key,value)
    lastopt = i+1

for acc in sys.argv[(lastopt+1):]:
    g = gtc.getGlycan(acc)
    outfile = acc + "." + imageWriter.format()
    imageWriter.writeImage(g,outfile)
