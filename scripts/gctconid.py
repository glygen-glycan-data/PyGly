#!/bin/env python2

import sys, time, traceback, hashlib, os, os.path, glob, csv
from collections import defaultdict

import findpygly
from pygly.alignment import GlycanEqual, GlycanImageEqual
from pygly.GlycanFormatter import GlycoCTFormat, WURCS20Format, GlycanParseError
from pygly.GlycanResource import GlyTouCanNoPrefetch

gtc = GlyTouCanNoPrefetch(verbose=False)

wp = WURCS20Format()
gp = GlycoCTFormat()
glyeq = GlycanEqual()
glyimgeq = GlycanImageEqual()

idmapfilename = "residmap.txt"
idmaps = defaultdict(dict)
if os.path.exists(idmapfilename):
    for r in csv.DictReader(open(idmapfilename),dialect='excel-tab'):
        idmaps[r['Accession']][int(r['GlycoCTResidueIndex'])] = r['CanonicalResidueIndex']

for acc in sys.argv[1:]:

    if acc.endswith('.txt'):
        gct = open(acc).read()
        acc = acc.split('.',1)[0]
    else:
        gct = gtc.getseq(acc,'glycoct')
        if not gct:
            gct = gtc.glycoct(acc)
    wcs = gtc.getseq(acc,'wurcs')

    if not wcs or not gct:
        continue

    try:
        gctgly = gp.toGlycan(gct)
        wcsgly = wp.toGlycan(wcs)
    except GlycanParseError:
        continue
    except:
        traceback.print_exc()
        print >>sys.stderr, acc
        sys.exit(1)

    idmap = []
    if glyeq.eq(gctgly,wcsgly,idmap=idmap):
        for gctmono,wcsmono in idmap:
            for gctid,wcsid in glyeq.monoidmap(gctmono,wcsmono):
                idmaps[acc][int(gctid)] = wcsid
    elif glyimgeq.eq(gctgly,wcsgly,idmap=idmap):
        for gctmono,wcsmono in idmap:
            for gctid,wcsid in glyimgeq.monoidmap(gctmono,wcsmono):
                idmaps[acc][int(gctid)] = wcsid

    filename = acc + ".txt"

    print >>sys.stderr, acc
    wh = open(filename,'w')
    wh.write(gct)
    wh.close()

wh = open(idmapfilename,'w')
print >>wh, "\t".join(map(str,["Accession","GlycoCTResidueIndex","CanonicalResidueIndex"]))
for acc in sorted(idmaps):
    for gctid in sorted(idmaps[acc]):
        print >>wh, "\t".join(map(str,[acc,gctid,idmaps[acc][gctid]]))
wh.close()
