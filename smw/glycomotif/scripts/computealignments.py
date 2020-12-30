#!/bin/env python27
import os
import sys
import time
import findpygly
import pygly.alignment
from pygly.GlycanFormatter import GlycoCTFormat, WURCS20Format
from pygly.GlycanResource.GlyTouCan import GlyTouCanNoCache

from getwiki import GlycoMotifWiki
w = GlycoMotifWiki()

if len(sys.argv) > 1:
    res_file_path = sys.argv[1] # "../data/motif_alignments.tsv"
else:
    res_file_path = None

wp = WURCS20Format()
gp = GlycoCTFormat()
gtc = GlyTouCanNoCache()

nodes_cache = pygly.alignment.ConnectedNodesCache()
loose_matcher = pygly.alignment.MotifInclusive(connected_nodes_cache=nodes_cache)
normal_matcher = pygly.alignment.MotifStrict(connected_nodes_cache=nodes_cache)
nred_matcher = pygly.alignment.NonReducingEndMotifStrict(connected_nodes_cache=nodes_cache)


motif_gobjs = {}
for m in w.itermotif():

    acc = m.get("glytoucan")
    if acc in motif_gobjs:
        continue

    try:
        motif_gobjs[acc] = wp.toGlycan(str(m.get("wurcs")))
    except:

        try:
            motif_gobjs[acc] = gp.toGlycan(m.get("glycoct"))
        except:
            continue


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

for glycan_acc, f, s in gtc.allseq(format="wurcs"):
    i += 1
    per = i / l

    nodes_cache.clear()
    try:
        glycan_obj = wp.toGlycan(s)
    except:
        continue

    if per > lastper:
        lastper += 0.1
        lapsed = time.time() - start_ts
        print "%0.2f Percent finished after %s, estimate %s remaining" % (per, secondtostr(lapsed), secondtostr(lapsed/per*(100-per)))


    for motif_acc in motif_accs:
        motif_gobj = motif_gobjs[motif_acc]

        loose = loose_matcher.leq(motif_gobj, glycan_obj, rootOnly=False, anywhereExceptRoot=False, strictMatch=False)

        if not loose:
            continue

        core = normal_matcher.leq(motif_gobj, glycan_obj, rootOnly=True, anywhereExceptRoot=False, strictMatch=True)
        substructure_partial = None
        if not core:
            substructure_partial = normal_matcher.leq(motif_gobj, glycan_obj, rootOnly=False, anywhereExceptRoot=True, strictMatch=True)

        whole = False
        if core and (len(list(motif_gobj.all_nodes())) == len(list(glycan_obj.all_nodes()))):
            whole = True

        nred = False
        if (core or substructure_partial):
            nred = nred_matcher.leq(motif_gobj, glycan_obj, strictMatch=True)

        substructure = core or substructure_partial

        res0 = [loose, core, substructure, whole, nred]

        if True not in res0:
            continue

        res0 = map(lambda x: "Y" if x else "N", res0)

        line = "\t".join([motif_acc, glycan_acc] + res0)

        # f1.write(line + "\n")
        result.append(line)

result = sorted(result)
if res_file_path:
    result_file = open(res_file_path, "w")
else:
    result_file = sys.stdout
result_file.write("Motif\tStructure\tLoose\tCore\tSubstructure\tWhole\tNon_Red\n")
result_file.write("\n".join(result))
if res_file_path:
    result_file.close()















