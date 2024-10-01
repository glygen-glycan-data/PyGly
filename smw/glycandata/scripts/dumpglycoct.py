#!/bin/env python3.12

import sys, time, traceback, hashlib, os, os.path, glob, csv
from collections import defaultdict

import findpygly
from pygly.alignment import GlycanEqual, GlycanImageEqual
from pygly.GlycanFormatter import GlycoCTFormat, WURCS20Format, GlycanParseError

from getwiki import GlycanData
w = GlycanData()

wp = WURCS20Format()
gp = GlycoCTFormat()
glyeq = GlycanEqual()
glyimgeq = GlycanImageEqual()

allfn = set(glob.glob(sys.argv[1]+"/G*.txt"))

idmapfilename = sys.argv[1] + "/residmap.txt"
idmaps = defaultdict(dict)
if os.path.exists(idmapfilename):
    for r in csv.DictReader(open(idmapfilename),dialect='excel-tab'):
        idmaps[r['Accession']][r['GlycoCTResidueIndex']] = r['CanonicalResidueIndex']

def check_idmap(gly1,gly2,idmap):
    ids1 = gly1.external_descriptor_ids()
    ids2 = gly2.external_descriptor_ids()
    idm = list(idmap.items())
    idm1 = [ p[0] for p in idm ]
    idm2 = [ p[1] for p in idm ]
    assert len(ids1) == len(set(ids1))
    assert len(ids2) == len(set(ids2))
    assert len(idm1) == len(set(idm1))
    assert len(idm2) == len(set(idm2))
    assert set(ids1) == set(idm1)
    assert set(ids2) == set(idm2)

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
                if value in ("GGM.001006","GGM.001010","GGM.001014","GGM.001016","GGM.001018","GGM.001033"):
                    inclass = True
                    break
                if value in ("GGM.001034",):
                    try:
                        monocnt = int(g.get_annotation_value(property="MonosaccharideCount",type="MonosaccharideCount",source="EdwardsLab"))
                    except ValueError:
                        monocnt = 0
                    try:
                        neuaccnt = int(g.get_annotation_value(property="NeuAcCount",type="MonosaccharideCount",source="EdwardsLab"))
                    except LookupError:
                        neuaccnt = 0
                    if monocnt == 1:
                        inclass = True
                        break
                    if monocnt == 2 and neuaccnt == 1:
                        inclass = True
                        break
    else:
        raise RuntimeError("Bad glycan-class...")

    if gct == None or wcs == None or not inclass:
        continue

    try:
        gctgly = gp.toGlycan(gct)
        wcsgly = wp.toGlycan(wcs)
    except GlycanParseError:
        continue
    except:
        traceback.print_exc()
        print(acc)
        sys.exit(1)

    idmap = []
    if glyeq.eq(gctgly,wcsgly,idmap=idmap):
        for gctmono,wcsmono in idmap:
            for gctid,wcsid in glyeq.monoidmap(gctmono,wcsmono):
                idmaps[acc][str(gctid)] = str(wcsid)
    elif glyimgeq.eq(gctgly,wcsgly,idmap=idmap):
        for gctmono,wcsmono in idmap:
            for gctid,wcsid in glyimgeq.monoidmap(gctmono,wcsmono):
                idmaps[acc][str(gctid)] = str(wcsid)

    check_idmap(gctgly,wcsgly,idmaps[acc])

    filename = sys.argv[1] + "/" + acc + ".txt"
    print(acc,file=sys.stderr)

    if filename in allfn:
        allfn.remove(filename)
        continue

    wh = open(filename,'w')
    wh.write(gct)
    wh.close()

for fn in allfn:
    acc = os.path.split(fn)[1][:-4]
    print("Removing:",acc,file=sys.stderr)
    os.unlink(fn)
    del idmaps[acc]

wh = open(idmapfilename,'w')
print("\t".join(map(str,["Accession","GlycoCTResidueIndex","CanonicalResidueIndex"])),file=wh)
for acc in sorted(idmaps):
    for gctid in sorted(idmaps[acc]):
        print("\t".join(map(str,[acc,gctid,idmaps[acc][gctid]])),file=wh)
wh.close()
