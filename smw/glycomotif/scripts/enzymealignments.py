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
from pygly.GlycanResource import GlyTouCan, GlyTouCanNoCache, GlycoTreeSandbox, GlycoTreeSandboxDev
from operator import itemgetter

from getwiki import GlycoMotifWiki

w = GlycoMotifWiki()

gtc = GlyTouCanNoCache()
# gts = GlycoTreeSandbox()
gts = GlycoTreeSandboxDev()

nodes_cache = pygly.alignment.ConnectedNodesCache()
strict_matcher = pygly.alignment.MotifStrict(connected_nodes_cache=nodes_cache)
strict_nred_matcher = pygly.alignment.NonReducingEndMotifStrict(connected_nodes_cache=nodes_cache)

accs = sys.argv[1:]

def motifgtc(accs):
    if len(accs) == 0:
       seen = set()
       for m in w.itermotif():
           gtcacc = m.get('glytoucan')
           if gtcacc in seen:
               continue
           seen.add(gtcacc)
           yield gtcacc
    else:
       for acc in accs:
           yield acc

motifs = {} 
motifids = {}
for macc in motifgtc(accs):
    # print >>sys.stderr, "Motif: ", macc
    motifs[macc] = gtc.getGlycan(macc,format='wurcs')
    if not motifs[macc] or motifs[macc].repeated():
        del motifs[macc]
        continue
    motifids[macc] = motifs[macc].external_descriptor_ids()
    assert len(motifids[macc]) == len(set(motifids[macc]))
    motifids[macc] = set(motifids[macc])

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
    for idmapids in idmaps:
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

def idmaps_toids(idmaps,aligner):
    newidmaps = []
    for idmap in idmaps:
        idmapids = [ ti for t in idmap for ti in aligner.monoidmap(*t) ]
        newidmaps.append(idmapids)
    return newidmaps

def idmaps_tosupportids(idmaps):
    retval = []
    for idmap in idmaps:
        supportids = []
        for m in support(idmap):
            supportids.extend(m.external_descriptor_ids())
        retval.append(supportids)
    return retval

print("\t".join(["Motif","MotifResidue","Structure","StructureResidue","AlignmentType","AlignmentIndex","GlycoTreeID","HumanEnzyme","MouseEnzyme"]))

def printres(macc,gacc,idmaps,supports,typ):
    for ind,(idmap,support) in enumerate(zip(idmaps,supports)):
        for mmid,gmid in idmap:
            if not sandbox.get(gmid):
                continue
            print "\t".join(map(str,[macc,mmid,gacc,gmid,typ,ind+1,
                                     sandbox.get(gmid,{}).get('resid','-'),
                                     sandbox.get(gmid,{}).get('humenz','-'),
                                     sandbox.get(gmid,{}).get('musenz','-'),
                                    ]))
        for gmid in support:
            if not sandbox.get(gmid):
                continue
            print "\t".join(map(str,[macc,"-",gacc,gmid,typ,ind+1,
                                     sandbox.get(gmid,{}).get('resid','-'),
                                     sandbox.get(gmid,{}).get('humenz','-'),
                                     sandbox.get(gmid,{}).get('musenz','-'),
                                    ]))

for gtsselector in ['mapped_N','mapped_O']:
  for enzdata in gts.allglycans(gtsselector):
    gacc = enzdata['accession']
    glycan = gtc.getGlycan(gacc,format='wurcs')
    if glycan.repeated():
        continue
    glycanids = glycan.external_descriptor_ids()
    assert len(glycanids) == len(set(glycanids))
    glycanids = set(glycanids)
    # print >>sys.stderr, "Structure:", gacc
    sandbox = dict()
    if len(enzdata.get('caveats',[])) > 0:
        continue
    for r in enzdata["residues"]:
        if r.get('glycotree',"none") == "none":
            continue
        gctind = r['canonical_residue_index']
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
        idmaps_core = []; idmaps_withids_core = []; support_withids_core = [];
        strict_core = strict_matcher.leq(motif, glycan, rootOnly=True, anywhereExceptRoot=False, underterminedLinkage=False, idmaps=idmaps_core)
        if strict_core:
            idmaps_withids_core = idmaps_toids(idmaps_core,strict_matcher)
            support_withids_core = idmaps_tosupportids(idmaps_core)
            if not check_idmaps(motifids[macc],glycanids,idmaps_withids_core):
                print >>sys.stderr, "Glycan:",gacc, "Motif:", macc
                sys.exit(1)

        idmaps_noncore = []; idmaps_withids_noncore = []; support_withids_noncore = [];
        strict_substructure_partial = strict_matcher.leq(motif, glycan, rootOnly=False, anywhereExceptRoot=True, underterminedLinkage=False, idmaps=idmaps_noncore)

        if strict_substructure_partial:
            idmaps_withids_noncore = idmaps_toids(idmaps_noncore,strict_matcher)
            support_withids_noncore = idmaps_tosupportids(idmaps_noncore)
            if not check_idmaps(motifids[macc],glycanids,idmaps_withids_noncore):
                print >>sys.stderr, "Glycan:",gacc, "Motif:", macc
                sys.exit(1)

        idmaps_whole = []; idmaps_withids_whole = []; support_withids_whole = [];
        if strict_core and strict_matcher.whole_glycan_match_check(motif, glycan):
            strict_whole = True
            idmaps_whole = list(idmaps_core)
            idmaps_withids_whole = list(idmaps_withids_core)
            support_withids_whole = list(support_withids_core)

        strict_substructure = (strict_core or strict_substructure_partial)
        idmaps_substructure = (idmaps_core + idmaps_noncore)
        idmaps_withids_substructure = (idmaps_withids_core + idmaps_withids_noncore)
        support_withids_substructure = (support_withids_core + support_withids_noncore)

        if strict_substructure and not check_idmaps(motifids[macc],glycanids,idmaps_withids_substructure):
            print >>sys.stderr, "Glycan:",gacc, "Motif:", macc
            sys.exit(1)

        idmaps_nonred = []; idmaps_withids_nonred = []; support_withids_nonred = [];
        if strict_substructure:
            strict_nred = strict_nred_matcher.leq(motif, glycan, underterminedLinkage=False, idmaps=idmaps_nonred)
            if strict_nred:
                idmaps_withids_nonred = idmaps_toids(idmaps_nonred,strict_matcher)
                support_withids_nonred = idmaps_tosupportids(idmaps_nonred)
                if not check_idmaps(motifids[macc],glycanids,idmaps_withids_nonred):
                    print >>sys.stderr, "Glycan:",gacc, "Motif:", macc
                    sys.exit(1)
        
        if strict_core:
            printres(macc,gacc,idmaps_withids_core,support_withids_core,"core")
        if strict_substructure:
            printres(macc,gacc,idmaps_withids_substructure,support_withids_substructure,"substr")
        if strict_nred:
            printres(macc,gacc,idmaps_withids_nonred,support_withids_nonred,"nonred")
        if strict_whole:
            printres(macc,gacc,idmaps_withids_whole,support_withids_whole,"whole")
