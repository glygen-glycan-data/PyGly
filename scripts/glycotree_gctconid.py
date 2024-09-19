#!/bin/env python3

import sys, time, traceback, hashlib, os, os.path, glob, csv
from collections import defaultdict

import findpygly
from pygly.alignment import GlycanEqual, GlycanImageEqual
from pygly.Monosaccharide import Mod
from pygly.GlycanFormatter import GlycoCTFormat, WURCS20Format, GlycanParseError
from pygly.GlycanResource import GlyTouCanNoCache, GlyTouCanNoPrefetch
from pygly.GlycanResource import GlycoMotif, GlycoMotifDev, GlycoMotifNoCache, GlycoMotifDevNoCache

gtc = GlyTouCanNoCache(verbose=False)
gm = GlycoMotifNoCache(local=True,verbose=False)

wp = WURCS20Format()
gp = GlycoCTFormat()
glyeq = GlycanEqual()
glyimgeq = GlycanImageEqual()

glycanclass = sys.argv[1]
assert glycanclass in ("N-linked","O-linked")
outdir = sys.argv[2]

if not os.path.exists(outdir):
    os.makedirs(outdir)

includeonly = None
if len(sys.argv) > 3:
    includeonly=set(sys.argv[3:])

allfn = set(glob.glob(outdir+"/G*.txt"))

idmapfilename = outdir + "/residmap.txt"
idmaps = defaultdict(dict)
if os.path.exists(idmapfilename):
    for r in csv.DictReader(open(idmapfilename),dialect='excel-tab'):
        idmaps[r['Accession']][r['GlycoCTResidueIndex']] = r['CanonicalResidueIndex']

validresidues = set(filter(None,"""
GlcNAc
Glc
GlcA
Man
Gal
GalNAc
Fuc
Xyl
Kdn
NeuAc
NeuGc
P
S
aldi
Count
""".split()))

def validcomp(comp,root):
    for k,v in comp.items():
        if k not in validresidues and v > 0:
            return False
    if comp['aldi'] > 1:
        return False
    if root.count_mod(Mod.aldi) == 0 and comp['aldi'] > 0:
        return False
    return True

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

def iternlinkedaccs():
    seen = set()
    for acc,strict,resids,linkids in gm.getstruct('GGM','001001'):
        if acc in seen:
            continue
        if includeonly and acc not in includeonly:
            continue
        gly = gtc.getGlycan(acc)
        if not gly or gly.repeated():
            continue
        comp = gly.iupac_composition(aggregate_basecomposition=False)
        # print(comp)
        if not validcomp(comp,gly.root()):
            continue
        yield acc
        seen.add(acc)

def iterolinkedaccs():
    seen = set()
    olinkedcores = [ "001006", "001010",
                     "001014", "001016", 
                     "001018", "001033",
                     "001034", "001035",
                     "001036", "001028",
                     "001031"
                   ]
    for olc in olinkedcores:
        for acc,strict,resids,linkids in gm.getstruct('GGM',olc):
            if acc in seen:
                continue
            if includeonly and acc not in includeonly:
                continue
            gly = gtc.getGlycan(acc)
            if not gly or gly.repeated():
                continue
            comp = gly.iupac_composition(aggregate_basecomposition=False)
            # print(acc,comp)
            if not validcomp(comp,gly.root()):
                continue
            if olc in "001034":
                if comp['Count'] not in (1,2):
                    continue
                if comp['Count'] == 2 and comp['NeuAc'] != 1:
                    continue
            yield acc
            seen.add(acc)

def iteracc():
    if glycanclass == "N-linked":
        return iternlinkedaccs()
    elif glycanclass == "O-linked":
        return iterolinkedaccs()
    return []

for acc in iteracc():
    
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
        print(acc)
        sys.exit(1)

    if wcsgly.repeated() or gctgly.repeated():
        continue

    if not wcsgly.has_root() or not gctgly.has_root():
        continue

    idmap = []
    if glyeq.eq(gctgly,wcsgly,idmap=idmap):
        for gctmono,wcsmono in idmap:
            for gctid,wcsid in glyeq.monoidmap(gctmono,wcsmono):
                idmaps[acc][str(gctid)] = str(wcsid)
        # print(acc,idmaps[acc])
    elif glyimgeq.eq(gctgly,wcsgly,idmap=idmap):
        for gctmono,wcsmono in idmap:
            for gctid,wcsid in glyimgeq.monoidmap(gctmono,wcsmono):
                idmaps[acc][str(gctid)] = str(wcsid)
        # print(acc,idmaps[acc])
    else:
        continue

    check_idmap(gctgly,wcsgly,idmaps[acc])

    filename = outdir + "/" + acc + ".txt"
    print(acc,file=sys.stderr)

    if filename in allfn:
        allfn.remove(filename)
        continue

    # print >>sys.stderr, acc
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
