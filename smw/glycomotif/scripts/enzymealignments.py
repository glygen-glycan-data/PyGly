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

wp = WURCS20Format()
gp = GlycoCTFormat()

nodes_cache = pygly.alignment.ConnectedNodesCache()

strict_matcher = pygly.alignment.MotifStrict(connected_nodes_cache=nodes_cache)
strict_nred_matcher = pygly.alignment.NonReducingEndMotifStrict(connected_nodes_cache=nodes_cache)

motifs = {}                                                                                                             
for mfn in sorted(glob.glob(sys.argv[1]+"/*.txt")):
    macc = mfn.rsplit('/',1)[-1][:-4]
    try:
        motifs[macc] = wp.toGlycan(open(mfn).read())
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
        strict_core, strict_substructure_partial, strict_whole, strict_nred = False, False, False, False
        idmaps_core = []
        strict_core = strict_matcher.leq(motif, glycan, rootOnly=True, anywhereExceptRoot=False, underterminedLinkage=False, idmaps=idmaps_core)
        idmaps_noncore = []
        strict_substructure_partial = strict_matcher.leq(motif, glycan, rootOnly=False, anywhereExceptRoot=True, underterminedLinkage=False, idmaps=idmaps_noncore)

        if strict_core and strict_matcher.whole_glycan_match_check(motif, glycan):
            strict_whole = True
            idmaps_whole = list(idmaps_core)

        strict_substructure = (strict_core or strict_substructure_partial)
        idmaps_substructure = list(idmaps_core) + list(idmaps_noncore)

        idmaps_nonred = []
        if strict_substructure:
            strict_nred = strict_nred_matcher.leq(motif, glycan, underterminedLinkage=False, idmaps=idmaps_nonred)
        
        if strict_core:
            printres(macc,gacc,idmaps_core,"core")
        if strict_substructure:
            printres(macc,gacc,idmaps_substructure,"substr")
        if strict_nred:
            printres(macc,gacc,idmaps_nonred,"nonred")
        if strict_whole:
            printres(macc,gacc,idmaps_whole,"whole")
