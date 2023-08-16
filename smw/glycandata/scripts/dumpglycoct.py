#!/bin/env python2

import sys, time, traceback, hashlib, os, os.path, glob, csv
from collections import defaultdict

import findpygly
from pygly.alignment import GlycanEqual
from pygly.GlycanFormatter import GlycoCTFormat, WURCS20Format, GlycanParseError

from getwiki import GlycanData
w = GlycanData()

wp = WURCS20Format()
gp = GlycoCTFormat()
glyeq = GlycanEqual()

allfn = set(glob.glob(sys.argv[1]+"/G*.txt"))

idmapfilename = sys.argv[1] + "/residmap.txt"
idmaps = defaultdict(dict)
if os.path.exists(idmapfilename):
    for r in csv.DictReader(open(idmapfilename),dialect='excel-tab'):
        idmaps[r['Accession']][int(r['GlycoCTResidueIndex'])] = int(r['CanonicalResidueIndex'])

glycanclass = sys.argv[2]
assert glycanclass in ("N-linked","O-linked")

for g in w.iterglycan():
    acc = g.get('accession')

    gct = None
    if g.has_annotations(property='GlycoCT',type='Sequence'):
	gct = g.get_annotation_value(property='GlycoCT',type='Sequence')

    wcs = None
    if g.has_annotations(property='WURCS',type='Sequence'):
	wcs = g.get_annotation_value(property='WURCS',type='Sequence')

    inclass = False
    if glycanclass == "N-linked":
        if g.has_annotations(property='ClassMotif',type='Motif',source='GlycoMotif'):
	    for value in g.get_annotation_values(property='ClassMotif',type='Motif',source='GlycoMotif'):
                if value == "GGM.001001":
		    inclass = True
		    break

    elif glycanclass == "O-linked":
        if g.has_annotations(property='ClassMotif',type='Motif',source='GlycoMotif'):
	    for value in g.get_annotation_values(property='ClassMotif',type='Motif',source='GlycoMotif'):
                if value in ("GGM.001006","GGM.001010","GGM.001014","GGM.001016","GGM.001018"):
		    inclass = True
		    break
    else:
	raise RuntimeError("Bad glycan-class...")

    if gct == None or not inclass:
	continue

    try:
        gctgly = gp.toGlycan(gct)
        wcsgly = wp.toGlycan(wcs)
    except GlycanParseError:
        continue
    except:
        traceback.print_exc()
        print acc
        sys.exit(1)

    idmap = []
    if not glyeq.eq(gctgly,wcsgly,idmap=idmap):
	continue

    for gctmono,wcsmono in idmap:
        idmaps[acc][int(gctmono.id())] = int(wcsmono.external_descriptor_id())

    filename = sys.argv[1] + "/" + acc + ".txt"

    if filename in allfn:
        allfn.remove(filename)	
	continue

    print >>sys.stderr, acc
    wh = open(filename,'w')
    wh.write(gct)
    wh.close()

for fn in allfn:
    acc = os.path.split(fn)[1][:-4]
    print >>sys.stderr, "Removing:",acc
    os.unlink(fn)
    del idmaps[acc]

wh = open(idmapfilename,'w')
print >>wh, "\t".join(map(str,["Accession","GlycoCTResidueIndex","CanonicalResidueIndex"]))
for acc in sorted(idmaps):
    for gctid in sorted(idmaps[acc]):
        print >>wh, "\t".join(map(str,[acc,gctid,idmaps[acc][gctid]]))
wh.close()
