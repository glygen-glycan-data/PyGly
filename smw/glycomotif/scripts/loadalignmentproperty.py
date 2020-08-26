#!/bin/env python27

import sys
from getwiki import GlycoMotifWiki
w = GlycoMotifWiki()


for m in w.itermotif():
    collection = m.get("collection")
    pageid = m.get("id")
    acc = m.get("glytoucan")
    redend = m.get("redend")

    if redend == None:
        redend = []

    # "Substructure" "Core" "Whole" "Nonreducing-End"
    alignment_type = []
    if True in redend:
        alignment_type.append("Core")

    if False in redend:
        alignment_type.append("Substructure")

    if redend == None or len(redend) == 0:
        alignment_type = ["Core", "Substructure"]

    if collection in ["GGM"]:

        try:
            idnum = pageid.split(".")[1]
            idnum = int(idnum)
            if idnum >= 60 and idnum <= 122:
                alignment_type = ["Whole-Glycan"]
        except:
            pass


        if pageid == "GGM.000034":
            alignment_type = ["Nonreducing-End"]

    print >> sys.stderr, "%s set to %s" % (m.get("id"), alignment_type)
    m.set("alignment", alignment_type)
    w.update(m)

