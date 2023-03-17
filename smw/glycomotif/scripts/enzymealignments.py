#!/bin/env python2
import os
import sys
import glob
import time
import json
from collections import defaultdict
import findpygly
import pygly.alignment
from pygly.GlycanFormatter import GlycoCTFormat, WURCS20Format, GlycanParseError
from operator import itemgetter

wp = WURCS20Format()
gp = GlycoCTFormat()

nodes_cache = pygly.alignment.ConnectedNodesCache()

strict_matcher = pygly.alignment.MotifStrict(connected_nodes_cache=nodes_cache)
strict_nred_matcher = pygly.alignment.NonReducingEndMotifStrict(connected_nodes_cache=nodes_cache)

motifs = {}                                                                                                             
motifids = {}
for mfn in sorted(glob.glob(sys.argv[1]+"/*.txt")):
    macc = mfn.rsplit('/',1)[-1][:-4]
    try:
        motifs[macc] = wp.toGlycan(open(mfn).read())
        motifids[macc] = [ m.id() for m in motifs[macc].all_nodes() ]
        assert len(motifids[macc]) == len(set(motifids[macc]))
        motifids[macc] = set(motifids[macc])
    except GlycanParseError:
        pass


def support(idmap):
    s = set()
    for mm,gm in idmap:
        while True:
            pl = gm.any_parent_link()
            if not pl:
                break
            gm = pl.parent()
            s.add(gm)
    for mm,gm in idmap:
        if gm in s:
            s.remove(gm)
    return s

def check_idmaps(motifids,glycanids,idmaps):
    bad = 0
    for idmap in idmaps:
        idmapids = [ (t[0].id(),t[1].id()) for t in idmap ]
        bad = 0
        if set(map(itemgetter(0),idmapids)) != motifids:
            bad += 1
        if not set(map(itemgetter(1),idmapids)) <= glycanids:
            print >>sys.stderr, set(map(itemgetter(1),idmapids))
            print >>sys.stderr, glycanids
            print >>sys.stderr, set(map(itemgetter(1),idmapids)) <= glycanids
            bad += 2
        if len(set(map(itemgetter(1),idmapids))) != len(motifids):
            bad += 4
        if len(map(itemgetter(0),idmapids)) != len(set(map(itemgetter(0),idmapids))):
            bad += 8
        if len(map(itemgetter(1),idmapids)) != len(set(map(itemgetter(1),idmapids))):
            bad += 16
        if len(set(map(itemgetter(0),idmapids))) != len(set(idmapids)):
            bad += 32
        if bad > 0:
            print >>sys.stderr, "Bad idmap:",bad,idmapids
            break
    return bad == 0

print("\t".join(["Motif","MotifResidue","Structure","StructureResidue","AlignmentType","AlignmentIndex","GlycoTreeID","HumanEnzyme","MouseEnzyme"]))

def printres(macc,gacc,idmaps,typ):
    for ind,idmap in enumerate(idmaps):
        for mm,gm in idmap:
            print "\t".join(map(str,[macc,mm.id(),gacc,gm.id(),typ,ind+1,
                                     sandbox.get(gm.id(),{}).get('resid','-'),
                                     sandbox.get(gm.id(),{}).get('humenz','-'),
                                     sandbox.get(gm.id(),{}).get('musenz','-'),
                                    ]))
        for gm in support(idmap):
            print "\t".join(map(str,[macc,"-",gacc,gm.id(),typ,ind+1,
                                     sandbox.get(gm.id(),{}).get('resid','-'),
                                     sandbox.get(gm.id(),{}).get('humenz','-'),
                                     sandbox.get(gm.id(),{}).get('musenz','-'),
                                    ]))

for gsfn in sorted(glob.glob(sys.argv[2]+"/*.txt")):
    gacc = gsfn.rsplit('/',1)[-1][:-4]
    glycan = gp.toGlycan(open(gsfn).read())
    glycanids = [ m.id() for m in glycan.all_nodes() ]
    assert len(glycanids) == len(set(glycanids))
    glycanids = set(glycanids)
    # print >>sys.stderr, gacc
    gsjfn = gsfn.rsplit('.',1)[0]+'.json'
    sandbox = dict()
    if os.path.exists(gsjfn):
        enzdata = json.loads(open(gsjfn).read())
        for r in enzdata["residues"]:
            if r.get('glycotree',"none") == "none":
                continue
            gctind = r['glycoct_index']
            resid = r['residue_id']
            enz = defaultdict(list)
            for e in r['enzymes']:
                enz[e['species']].append(e['gene_name'])
            sandbox[gctind] = dict(resid=resid,
                                   humenz=",".join(enz['Homo sapiens'] if len(enz['Homo sapiens']) > 0 else "-"),
                                   musenz=",".join(enz['Mus musculus'] if len(enz['Mus musculus']) > 0 else "-"))

    for macc,motif in motifs.items():
        # print >>sys.stderr, " ",macc
        strict_core, strict_substructure, strict_substructure_partial, strict_whole, strict_nred = False, False, False, False, False
        idmaps_core = []
        strict_core = strict_matcher.leq(motif, glycan, rootOnly=True, anywhereExceptRoot=False, underterminedLinkage=False, idmaps=idmaps_core)
        if not check_idmaps(motifids[macc],glycanids,idmaps_core):
            print >>sys.stderr, "Glycan:",gacc, "Motif:", macc
            sys.exit(1)

        idmaps_noncore = []
        strict_substructure_partial = strict_matcher.leq(motif, glycan, rootOnly=False, anywhereExceptRoot=True, underterminedLinkage=False, idmaps=idmaps_noncore)

        if not check_idmaps(motifids[macc],glycanids,idmaps_noncore):
            print >>sys.stderr, "Glycan:",gacc, "Motif:", macc
            sys.exit(1)

        if strict_core and strict_matcher.whole_glycan_match_check(motif, glycan):
            strict_whole = True
            idmaps_whole = list(idmaps_core)

        strict_substructure = (strict_core or strict_substructure_partial)
        idmaps_substructure = (idmaps_core + idmaps_noncore)

        if not check_idmaps(motifids[macc],glycanids,idmaps_substructure):
            print >>sys.stderr, "Glycan:",gacc, "Motif:", macc
            sys.exit(1)

        idmaps_nonred = []
        if strict_substructure:
            strict_nred = strict_nred_matcher.leq(motif, glycan, underterminedLinkage=False, idmaps=idmaps_nonred)
            if not check_idmaps(motifids[macc],glycanids,idmaps_nonred):
                print >>sys.stderr, "Glycan:",gacc, "Motif:", macc
                sys.exit(1)
        
        if strict_core:
            printres(macc,gacc,idmaps_core,"core")
        if strict_substructure:
            printres(macc,gacc,idmaps_substructure,"substr")
        if strict_nred:
            printres(macc,gacc,idmaps_nonred,"nonred")
        if strict_whole:
            printres(macc,gacc,idmaps_whole,"whole")
