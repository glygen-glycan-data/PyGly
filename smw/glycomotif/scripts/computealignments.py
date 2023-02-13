#!/bin/env python2
import os
import sys
import time
import findpygly
import pygly.alignment
from pygly.GlycanFormatter import GlycoCTFormat, WURCS20Format
from pygly.GlycanResource.GlyTouCan import GlyTouCanNoCache, GlyTouCan
from pygly.GlycanResource.GlyCosmos import GlyCosmosNoCache, GlyCosmos
from getwiki import GlycoMotifWiki
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
if len(sys.argv) > 2:
  for acc in sys.argv[2:]:
    gly = gtc.getGlycan(acc)
    if gly:
        motif_gobjs[acc] = gly
else:
  for m in w.itermotif():
    acc = m.get("glytoucan")

    if acc in motif_gobjs:
        continue


    gly = None
    if not gly:
        try:
            motif_gobjs[acc] = wp.toGlycan(str(m.get("wurcs")))
        except:
            pass
    if not gly:
        try:
            motif_gobjs[acc] = gp.toGlycan(m.get("glycoct"))
        except:
            pass

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

for glycan_acc, f, s in gtc.allseq(format="wurcs"):
    i += 1
    per = i / l

    if glycan_acc in archived:
        continue

    nodes_cache.clear()
    try:
        glycan_obj = wp.toGlycan(s)
    except:
        continue

    if per > lastper:
        lastper += 0.1
        lapsed = time.time() - start_ts
        print >> sys.stderr, "%0.2f Percent finished after %s, estimate %s remaining" % (per, secondtostr(lapsed), secondtostr(lapsed/per*(100-per)))


    for motif_acc in motif_accs:
        motif_gobj = motif_gobjs[motif_acc]

        # Loose match first
        loose_core = loose_matcher.leq(motif_gobj, glycan_obj, rootOnly=True, anywhereExceptRoot=False, underterminedLinkage=True)
        loose_substructure_partial = False
        if not loose_core:
            loose_substructure_partial = loose_matcher.leq(motif_gobj, glycan_obj, rootOnly=False, anywhereExceptRoot=True, underterminedLinkage=True)

        loose_substructure = (loose_core or loose_substructure_partial)

        loose_whole = False
        if loose_core and loose_matcher.whole_glycan_match_check(motif_gobj, glycan_obj):
            loose_whole = True

        loose_nred = False
        if not motif_gobj.repeated() and not glycan_obj.repeated() and loose_substructure:
            loose_nred = loose_nred_matcher.leq(motif_gobj, glycan_obj, underterminedLinkage=True)


        # if inclusive, then try to match strict
        strict_core, strict_substructure, strict_substructure_partial, strict_whole, strict_nred = False, False, False, False, False
        if loose_core:
            strict_core = strict_matcher.leq(motif_gobj, glycan_obj, rootOnly=True, anywhereExceptRoot=False, underterminedLinkage=False)

        if loose_substructure:
            if not strict_core:
                strict_substructure_partial = strict_matcher.leq(motif_gobj, glycan_obj, rootOnly=False, anywhereExceptRoot=True, underterminedLinkage=False)

        if strict_core and strict_matcher.whole_glycan_match_check(motif_gobj, glycan_obj):
            strict_whole = True

        strict_substructure = (strict_core or strict_substructure_partial)

        if loose_nred and strict_substructure:
            strict_nred = strict_nred_matcher.leq(motif_gobj, glycan_obj, underterminedLinkage=False)



        res0 = [loose_core, loose_substructure, loose_whole, loose_nred,
                strict_core, strict_substructure, strict_whole, strict_nred]

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
result_file.write("Motif\tStructure\tCore_Inclusive\tSubstructure_Inclusive\tWhole_Inclusive\tNon_Red_Inclusive\tCore_Strict\tSubstructure_Strict\tWhole_Strict\tNon_Red_Strict\n")
result_file.write("\n".join(result))
if res_file_path:
    result_file.close()















