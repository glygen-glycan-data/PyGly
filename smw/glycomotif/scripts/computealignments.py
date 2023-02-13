
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
    
print(res_file_path)

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

def check_idmaps(motifids,glycanids,idmaps): ###ID MAPS
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


def get_match_index (idmaps):
    
    allstructids = set()
    for idmap in idmaps: 
        
        strucids = [t[1].id() for t in idmap]
  
        allstructids.update(strucids)
    x = "Y:"+",".join(str(i) for i in sorted(allstructids))
    
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
gtc = GlyTouCan()

nodes_cache = pygly.alignment.ConnectedNodesCache()

loose_matcher = pygly.alignment.MotifInclusive(connected_nodes_cache=nodes_cache)
loose_nred_matcher = pygly.alignment.NonReducingEndMotifInclusive(connected_nodes_cache=nodes_cache)

strict_matcher = pygly.alignment.MotifStrict(connected_nodes_cache=nodes_cache)
strict_nred_matcher = pygly.alignment.NonReducingEndMotifStrict(connected_nodes_cache=nodes_cache)

motif_gobjs = {}
motifids = {}


if len(sys.argv) > 2:
  for acc in sys.argv[2:]:
    print(acc)
    gly = gtc.getGlycan(acc)
    if acc in motif_gobjs:
        continue
    if gly:
        motif_gobjs[acc] = gly
        motifids[acc] = [ m.id() for m in motif_gobjs[acc].all_nodes() ]
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
            motifids[acc] = [ m.id() for m in motif_gobjs[acc].all_nodes() ]
            assert len(motifids[acc]) == len(set(motifids[acc]))
            motifids[acc] = set(motifids[acc]) 
 

archived = set()
gco = GlyCosmos()
for acc in gco.archived():
    acc = acc["accession"]
    archived.add(acc)

def secondtostr(i):
    i = int(i)

    h = i / 3600
    m = (i - h * 3600) / 60

    h = str(h)
    if len(h) == 1:
        h = "0" + h

    m = str(m)
    if len(m) == 1:
        m = "0" + m

    return "%sh:%sm" % (h, m)

# f1 = open("tmp.txt", "w")
result = []
i, l, lastper = 0.0, len(list(gtc.allseq(format="wurcs"))) / 100.0, 0


start_ts = time.time()
motif_accs = motif_gobjs.keys()

archived = set()
gco = GlyCosmosNoCache()
for acc in gco.archived():
    acc = acc["accession"]
    archived.add(acc)
    #print(acc)


def secondtostr(i):
    i = int(i)

    h = i / 3600
    m = (i - h * 3600) / 60

    h = str(h)
    if len(h) == 1:
        h = "0" + h

    m = str(m)
    if len(m) == 1:
        m = "0" + m

    return "%sh:%sm" % (h, m)

# f1 = open("tmp.txt", "w")
result = []
i, l, lastper = 0.0, len(list(gtc.allseq(format="wurcs"))) / 100.0, 0


start_ts = time.time()
motif_accs = motif_gobjs.keys()

final_results = []

for glycan_acc, f, s in gtc.allseq(format="wurcs"):
    #print(glycan_acc)
    
        
    i += 1
    per = i / l

    if glycan_acc in archived:
        continue
 
    nodes_cache.clear()
    
    try:
        glycan_obj = wp.toGlycan(s)
        glycan = wp.toGlycan(s)
        #print(glycan_obj)
        glycanids = [ m.id() for m in glycan_obj.all_nodes() ]
        assert len(glycanids) == len(set(glycanids))
        glycanids = set(glycanids)
    except:
        continue

    if per > lastper:
        lastper += 0.1
        lapsed = time.time() - start_ts
        print >> sys.stderr, "%0.2f Percent finished after %s, estimate %s remaining" % (per, secondtostr(lapsed), secondtostr(lapsed/per*(100-per)))


    for macc,motif in motif_gobjs.items():
        
        
        idmaps_loose_core = []
        idmaps_loose_substructure = []
        idmaps_loose_noncore = []
        idmaps_loose_whole = [] 
        idmaps_loose_nred = [] 
        idmaps_strict_core = []
        idmaps_strict_noncore =[]
        idmaps_strict_substructure = []
        idmaps_strict_substructure_partial=[]
        idmaps_strict_whole = []
        idmaps_strict_nred = []
        
       
        loose_core = loose_matcher.leq(motif, glycan, rootOnly=True, anywhereExceptRoot=False, underterminedLinkage=True, idmaps=idmaps_loose_core)
        
        loose_substructure_noncore = False

        
    

      

       
        loose_substructure_noncore = loose_matcher.leq(motif, glycan, rootOnly=False, anywhereExceptRoot=True, underterminedLinkage=True, idmaps=idmaps_loose_noncore)

        
         
        
        loose_substructure = loose_core or loose_substructure_noncore
        
        
        idmaps_loose_substructure = list(idmaps_loose_core) + list(idmaps_loose_noncore)

        loose_whole = False
        
        if loose_core and loose_matcher.whole_glycan_match_check(motif, glycan):

            loose_whole=True
            idmaps_loose_whole = list(idmaps_loose_core)
            
            
        
        loose_nred=False
       
        if not motif.repeated() and not glycan.repeated() and loose_substructure:
            loose_nred = loose_nred_matcher.leq(motif, glycan, underterminedLinkage=True, idmaps=idmaps_loose_nred)

    
        strict_core, strict_substructure_partial, strict_whole, strict_nred = False, False, False, False
        
       
        
        if loose_core:
            
            strict_core = strict_matcher.leq(motif, glycan, rootOnly=True, anywhereExceptRoot=False, underterminedLinkage=False, idmaps=idmaps_strict_core)
            
      
                
        strict_noncore = strict_matcher.leq(motif, glycan, rootOnly=False, anywhereExceptRoot=True, underterminedLinkage=False, idmaps=idmaps_strict_noncore) #partial changed here
     

        if strict_core and strict_matcher.whole_glycan_match_check(motif, glycan):
            strict_whole = True
            idmaps_strict_whole = list(idmaps_strict_core)
        
        strict_substructure = strict_core or strict_noncore
        idmaps_strict_substructure = list(idmaps_strict_core) + list(idmaps_strict_noncore) #partial changed here

       
        if loose_nred and strict_substructure:
            strict_nred = strict_nred_matcher.leq(motif, glycan, underterminedLinkage=False, idmaps=idmaps_strict_nred)
            
        
        
        res0 = [loose_core, loose_substructure, loose_whole, loose_nred, 
                strict_core, strict_substructure,strict_whole, strict_nred]
        
        
            

        
        if True not in res0: 
            continue

        index_numbers = [idmaps_loose_core,idmaps_loose_substructure,idmaps_loose_whole,idmaps_loose_nred, idmaps_strict_core, idmaps_strict_substructure , idmaps_strict_whole, idmaps_strict_nred]
        
        res1 = []
        
        
        
        
        x=0
        for mt,idmaps in zip(res0,index_numbers):
        
            if mt:
              
                indices=get_match_index(idmaps)
                
                res1.append(indices)
            else:
                res1.append("N")
            
        
        #print >>sys.stderr,"Glycan:",glycan_acc, "Motif:", macc , " ".join(str(n) for n in res1)

        line = "\t".join([macc, glycan_acc] + res1)
        
        
        #uncomment here for line by line check 
        
       
        # values = line.split("\t")
        # my_indices=values[2:]
        # my_check = alignmentindexchecker.set_checker(my_indices)
        # if my_check == 0:
        #     print("good")
        # else:
        #     print("invalid:", line)
        #     #sys.exit()

         
        

        # f1.write(line + "\n")
        result.append(line)

result = sorted(result)

if res_file_path:
    result_file = open(res_file_path, "w")
else:
    result_file = sys.stdout
result_file.write("Motif\tStructure\tCore_Inclusive\tSubstructure_Inclusive\tWhole_Inclusive\tNon_Red_Inclusive\tCore_Strict\tSubstructure_Strict\tWhole_Strict\tNon_Red_Strict\n")
result_file.write("\n".join(result))
if res_file_path:
    result_file.close()





