import os
import sys
import time
import findpygly
import pygly.alignment
from pygly.GlycanFormatter import GlycoCTFormat, WURCS20Format
from pygly.GlycanResource.GlyTouCan import GlyTouCanNoCache

res_file_path = sys.argv[1] # "../data/motif_alignments.tsv"

wp = WURCS20Format()
gp = GlycoCTFormat()
gtc = GlyTouCanNoCache()

nodes_cache = pygly.alignment.ConnectedNodesCache()
normal_matcher = pygly.alignment.GlyGenMotif(connected_nodes_cache=nodes_cache)
nred_matcher = pygly.alignment.GlyGenMotifNonReducingEnd(connected_nodes_cache=nodes_cache)
gtcm_gen = pygly.alignment.GlyTouCanMotif(connected_nodes_cache=nodes_cache)


from getwiki import GlycoMotifWiki
w = GlycoMotifWiki()
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


start_ts = time.time()
motif_accs = motif_gobjs.keys()

# f1 = open("tmp.txt", "w")
result = []
i, l, lastper = 0.0, 90000 / 100.0, 0
for glycan_acc, f, s in gtc.allseq(format="wurcs"):
    nodes_cache.clear()
    try:
        glycan_obj = wp.toGlycan(s)
    except:
        continue

    i += 1
    per = i / l
    if per > lastper:
        lastper += 1
        print "%0.2f Percent finished after %is" % (per, time.time() - start_ts)

    for motif_acc in motif_accs:
        motif_gobj = motif_gobjs[motif_acc]

        core = normal_matcher.leq(motif_gobj, glycan_obj, rootOnly=True, anywhereExceptRoot=False)
        # Save as much runtime as possible
        substructure_partial = None
        if not core:
            substructure_partial = normal_matcher.leq(motif_gobj, glycan_obj, rootOnly=False, anywhereExceptRoot=True)

        whole = False
        if core and (len(list(motif_gobj.all_nodes())) == len(list(glycan_obj.all_nodes()))):
            whole = gtcm_gen.leq(motif_gobj, glycan_obj, rootOnly=True, anywhereExceptRoot=False)

        nred = False
        if (core or substructure_partial):
            nred = nred_matcher.leq(motif_gobj, glycan_obj)

        substructure = core or substructure_partial

        res0 = [core, substructure, whole, nred]

        if True not in res0:
            continue

        res0 = map(lambda x: "Y" if x else "N", res0)

        line = "\t".join([motif_acc, glycan_acc] + res0)

        # f1.write(line + "\n")
        result.append(line)

result = sorted(result)
result_file = open(res_file_path, "w")
result_file.write("Motif\tStructure\tCore\tSubstructure\tWhole\tNon_Red\n")
result_file.write("\n".join(result))
result_file.close()















