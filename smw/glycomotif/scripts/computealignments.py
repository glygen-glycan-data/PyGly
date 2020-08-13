import os
import sys
import time
import findpygly
from pygly.alignment import *
from pygly.GlycanFormatter import GlycoCTFormat, WURCS20Format
from pygly.GlycanResource.GlyTouCan import GlyTouCan, GlyTouCanNoCache

res_file_path = sys.argv[1] # "../data/motif_alignments.tsv"

wp = WURCS20Format()
gp = GlycoCTFormat()
gtc = GlyTouCan()

nodes_cache = ConnectedNodesCache()
gtcm = GlyTouCanMotif(connected_nodes_cache=nodes_cache)
gtcmnre = GlyTouCanMotifNonReducingEnd(connected_nodes_cache=nodes_cache)
mm_st = MotifAllowOptionalSub(connected_nodes_cache=nodes_cache)


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

        c = mm_st.leq(motif_gobj, glycan_obj, rootOnly=True, anywhereExceptRoot=False)
        d = mm_st.leq(motif_gobj, glycan_obj, rootOnly=False, anywhereExceptRoot=True)

        a = False
        if c:
            a = gtcm.leq(motif_gobj, glycan_obj, rootOnly=True, anywhereExceptRoot=False)

        b = False
        if d:
            b = gtcm.leq(motif_gobj, glycan_obj, rootOnly=False, anywhereExceptRoot=True)

        e = False
        if a or b:
            e = gtcmnre.leq(motif_gobj, glycan_obj, rootOnly=False, anywhereExceptRoot=False)

        f = False
        if a and len(list(glycan_obj.all_nodes())) == len(list(motif_gobj.all_nodes())):
            f = True

        res0 = [a, b, c, d, e, f]

        if True not in res0:
            continue

        res0 = map(lambda x: "Y" if x else "N", res0)

        line = "\t".join([motif_acc, glycan_acc] + res0)

        # f1.write(line + "\n")
        result.append(line)

result = sorted(result)
result_file = open(res_file_path, "w")
result_file.write("motif\tstructure\tgeneral_red\tgeneral_other\tst_red\tst_other\tnon_red_only\twhole\n")
result_file.write("\n".join(result))















