#!/bin/env python3

import os
import sys
import time
import glob
import findpygly
import pygly.alignment
from pygly.GlycanFormatter import GlycoCTFormat, WURCS20Format
from pygly.GlycanResource.GlyTouCan import GlyTouCanNoCache, GlyTouCan
from pygly.GlycanResource.GlyCosmos import GlyCosmosNoCache, GlyCosmos
from getwiki import GlycoMotifWiki
import alignmentindexchecker

import json
from collections import defaultdict
from pygly.GlycanFormatter import GlycoCTFormat, WURCS20Format, GlycanParseError
from operator import itemgetter


result = []

if len(sys.argv) > 1:
    res_file_path = sys.argv[1] # "../data/motif_alignments.tsv"
else:
    res_file_path = None

if res_file_path == "-":
    res_file_path = None
    
def check_idmaps(motifids,glycanids,idmaps): ###ID MAPS
    bad = 0
    for idmap in idmaps:
        idmapids = [ (t[0].id(),t[1].id()) for t in idmap ]
        bad = 0
        if set(map(itemgetter(0),idmapids)) != motifids:
            bad += 1
        if not set(map(itemgetter(1),idmapids)) <= glycanids:
            print(set(map(itemgetter(1),idmapids)),file-sys.stderr)
            print(glycanids,file-sys.stderr)
            print(set(map(itemgetter(1),idmapids)) <= glycanids,file-sys.stderr)
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
            print("Bad idmap:",bad,idmapids,file=sys.stderr)
            break
    return bad == 0

def idmaps_toids(idmaps,aligner):
    newidmaps = []
    for idmap in idmaps:
        idmapids = [ ti for t in idmap for ti in aligner.monoidmap(*t) ]
        newidmaps.append(idmapids)
    return newidmaps

def get_match_index (motifacc,structacc,idmaps,glycan=None):
    allstructids = set()
    allstructlinkids = set()
    for idmap in idmaps: 
        structids = set(t[1] for t in idmap)
        allstructids.update(structids)
        monoids = set(t1 for t1 in structids if '.' not in t1)
        linkids = set()
        if glycan:
            for l in glycan.all_links(uninstantiated=True):
                if l.parent().external_descriptor_id() in monoids and \
                   l.child().external_descriptor_id() in monoids:
                    if l.instantiated():
                        linkids.add((l.parent().external_descriptor_id(),l.child().external_descriptor_id()))
                    else:
                        linkids.add(('',l.child().external_descriptor_id()))
            if len(monoids) != (len(linkids) + 1):
                print("Warning: Bad linkids length motifacc: %s structacc: %s"%(motifacc,structacc),file=sys.stderr)
            allstructlinkids.update(linkids)
    x = "Y:"+",".join(str(i) for i in sorted(allstructids))
    if glycan:
        x += ":"+",".join("%s-%s"%p for p in sorted(allstructlinkids))
    return(x)


w = GlycoMotifWiki()

if len(sys.argv) > 1:
    res_file_path = sys.argv[1] # "../data/motif_alignments.tsv"
else:
    res_file_path = None

if res_file_path == "-":
    res_file_path = None

wp = WURCS20Format()
gp = GlycoCTFormat()
gtc = GlyTouCanNoCache()

nodes_cache = pygly.alignment.ConnectedNodesCache()

loose_matcher = pygly.alignment.MotifInclusive(connected_nodes_cache=nodes_cache)
loose_nred_matcher = pygly.alignment.NonReducingEndMotifInclusive(connected_nodes_cache=nodes_cache)

strict_matcher = pygly.alignment.MotifStrict(connected_nodes_cache=nodes_cache)
strict_nred_matcher = pygly.alignment.NonReducingEndMotifStrict(connected_nodes_cache=nodes_cache)

motif_gobjs = {}
motifids = {}


if len(sys.argv) > 2:
  for acc in sys.argv[2:]:
    # print(acc)
    gly = gtc.getGlycan(acc)
    if acc in motif_gobjs:
        continue
    if gly:
        motif_gobjs[acc] = gly
        motifids[acc] = motif_gobjs[acc].external_descriptor_ids()
        assert len(motifids[acc]) == len(set(motifids[acc]))
        motifids[acc] = set(motifids[acc]) 
   
else:
      for m in w.itermotif():
        acc = m.get("glytoucan")
        
        if acc in motif_gobjs:
            continue
        gly = gtc.getGlycan(acc)
        if gly:
            motif_gobjs[acc] = gly
            motifids[acc] = motif_gobjs[acc].external_descriptor_ids()
            assert len(motifids[acc]) == len(set(motifids[acc]))
            motifids[acc] = set(motifids[acc]) 
 

archived = set()
gco = GlyCosmosNoCache()
for acc in gco.archived():
    acc = acc["accession"]
    archived.add(acc)

def secondtostr(i):
    i = int(i)

    h = i // 3600
    m = (i % 3600) // 60

    h = str(h)
    if len(h) == 1:
        h = "0" + h

    m = str(m)
    if len(m) == 1:
        m = "0" + m

    return "%sh:%sm" % (h, m)

# f1 = open("tmp.txt", "w")
gtcallseq = sorted(gtc.allseq(format="wurcs"))
result = []
i, l, lastper = 0.0, len(gtcallseq)) / 100.0, 0

start_ts = time.time()
motif_accs = motif_gobjs.keys()

final_results = []

if res_file_path:
    result_file = open(res_file_path, "w")
else:
    result_file = sys.stdout
result_file.write("Motif\tStructure\tCore_Inclusive\tSubstructure_Inclusive\tWhole_Inclusive\tNon_Red_Inclusive\tCore_Strict\tSubstructure_Strict\tWhole_Strict\tNon_Red_Strict\n")

for glycan_acc, f, s in gtcallseq:
    #print(glycan_acc)

    i += 1
    per = i / l

    # if glycan_acc in archived:
    #     continue
 
    nodes_cache.clear()
    
    try:
        glycan = wp.toGlycan(s)
        glycanids = glycan.external_descriptor_ids()
        assert len(glycanids) == len(set(glycanids))
        glycanids = set(glycanids)
    except:
        continue

    if per > lastper:
        lastper += 0.1
        lapsed = time.time() - start_ts
        print("%0.2f Percent finished after %s, estimate %s remaining" % (per, secondtostr(lapsed), secondtostr(lapsed/per*(100-per))),file=sys.stderr)


    for macc,motif in sorted(motif_gobjs.items()):
        
       
        idmaps_loose_core = []
        loose_core = loose_matcher.leq(motif, glycan, rootOnly=True, anywhereExceptRoot=False, underterminedLinkage=True, idmaps=idmaps_loose_core)
        idmaps_loose_core = idmaps_toids(idmaps_loose_core,loose_matcher)
        
        idmaps_loose_noncore = []
        loose_noncore = loose_matcher.leq(motif, glycan, rootOnly=False, anywhereExceptRoot=True, underterminedLinkage=True, idmaps=idmaps_loose_noncore)
        idmaps_loose_noncore = idmaps_toids(idmaps_loose_noncore,loose_matcher)
 
        loose_substructure = loose_core or loose_noncore
        idmaps_loose_substructure = list(idmaps_loose_core) + list(idmaps_loose_noncore)

        loose_whole = False
        idmaps_loose_whole = [] 
        if loose_core and loose_matcher.whole_glycan_match_check(motif, glycan):
            loose_whole=True
            idmaps_loose_whole = list(idmaps_loose_core)
            
        loose_nred=False
        idmaps_loose_nred = [] 
        if not motif.repeated() and not glycan.repeated() and loose_substructure:
            loose_nred = loose_nred_matcher.leq(motif, glycan, underterminedLinkage=True, idmaps=idmaps_loose_nred)
            idmaps_loose_nred = idmaps_toids(idmaps_loose_nred,loose_nred_matcher)

        strict_core, strict_noncore, strict_whole, strict_nred = False, False, False, False
        idmaps_strict_core = []
        idmaps_strict_noncore =[]
        idmaps_strict_whole = []
        idmaps_strict_nred = []
        
        if loose_core:
            strict_core = strict_matcher.leq(motif, glycan, rootOnly=True, anywhereExceptRoot=False, underterminedLinkage=False, idmaps=idmaps_strict_core)
            idmaps_strict_core = idmaps_toids(idmaps_strict_core,strict_matcher)
      
        if loose_noncore:
            strict_noncore = strict_matcher.leq(motif, glycan, rootOnly=False, anywhereExceptRoot=True, underterminedLinkage=False, idmaps=idmaps_strict_noncore)
            idmaps_strict_noncore = idmaps_toids(idmaps_strict_noncore,strict_matcher)
     
        strict_substructure = strict_core or strict_noncore
        idmaps_strict_substructure = list(idmaps_strict_core) + list(idmaps_strict_noncore) #partial changed here

        if strict_core and strict_matcher.whole_glycan_match_check(motif, glycan):
            strict_whole = True
            idmaps_strict_whole = list(idmaps_strict_core)
       
        if loose_nred and strict_substructure:
            strict_nred = strict_nred_matcher.leq(motif, glycan, underterminedLinkage=False, idmaps=idmaps_strict_nred)
            idmaps_strict_nred = idmaps_toids(idmaps_strict_nred,strict_nred_matcher)
            
        res0 = [loose_core, loose_substructure, loose_whole, loose_nred, 
                strict_core, strict_substructure,strict_whole, strict_nred]
        
        if True not in res0: 
            continue

        index_numbers = [idmaps_loose_core,idmaps_loose_substructure,idmaps_loose_whole,idmaps_loose_nred, idmaps_strict_core, idmaps_strict_substructure , idmaps_strict_whole, idmaps_strict_nred]
        
        res1 = []
        
        for mt,idmaps in zip(res0,index_numbers):
            if mt:
                indices=get_match_index(macc,glycan_acc,idmaps,glycan)
                res1.append(indices)
            else:
                res1.append("N")
        
        #print >>sys.stderr,"Glycan:",glycan_acc, "Motif:", macc , " ".join(str(n) for n in res1)
        
        line = [macc, glycan_acc] + res1
        
        #uncomment here for line by line check 
        
       
        # my_check = alignmentindexchecker.set_checker(res1)
        # if my_check == 0:
        #     print("good")
        # else:
        #     print("invalid:", line)
        #     #sys.exit()

        print("\t".join(line),file=result_file)


if res_file_path:
    result_file.close()





