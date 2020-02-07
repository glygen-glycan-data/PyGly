#!/bin/env python27

from collections import defaultdict
import sys
import csv
import copy

motif_alignment_file = sys.argv[1]#"../data/motif_alignment.tsv"

motif_alignment = {}
motif_alignment_single = {
    "red_only": [],
    "other": []
}
with open(motif_alignment_file) as f:
    for r in csv.reader(f, delimiter="\t"):
        motif_acc, match, reducing_end_match = r

        if motif_acc not in motif_alignment:
            motif_alignment[motif_acc] = copy.deepcopy(motif_alignment_single)

        if reducing_end_match == "re":
            motif_alignment[motif_acc]["red_only"].append(match)

        if reducing_end_match == "other":
            motif_alignment[motif_acc]["other"].append(match)


from getwiki import GlycoMotifWiki, AllMotif
w = GlycoMotifWiki()


for m in w.itermotif():

    if m.get('collection') != AllMotif.id:
        continue

    motif_acc = m.get("glytoucan")

    #if motif_acc == "G00012MO":
    #    m.data["redend_alignments"] = ["test1"]
    #    m.data["other_alignments"] = ["testx"]

    #    if w.update(m):
    #        print >> sys.stderr, motif_acc

    #if motif_acc == "G00012MO":
    #    print m.data["redend_alignments"]
    #    print m.data["other_alignments"]

    if not motif_acc in motif_alignment:
        continue

    # TODO is there API for change data?
    m.data["redend_alignments"] = sorted(motif_alignment["red_only"])
    m.data["other_alignments"] = sorted(motif_alignment["other"])

    if w.update(m):
        print >> sys.stderr, motif_acc
