#!/bin/env python27

import os
import sys

from getwiki import GlycoMotifWiki
w = GlycoMotifWiki()


for m in w.itermotif():
    collection = m.get("collection")
    pageid = m.get("id")
    acc = m.get("glytoucan")
    redend = m.get("redend")

    idnum = pageid.split(",")[1]

    # CCRC CCRC.000114 G26199BG [True]
    # print collection, pageid, acc, redend

    if collection not in ["GGM"]:
        continue

    alignment_type = []
    # alignment_type.append("")
    # "Substructure" "Core" "Whole" "Non-Reducing"

    if True in redend:
        alignment_type.append("")
    else:
        alignment_type.append("")

    try:
        idnum = int(idnum)
        if idnum >=60 and idnum <=122:
            alignment_type = ["Whole"]
    except:
        pass


    if pageid == "GGM.000034":
        alignment_type = ["Non-Reducing"]


    print "%s set to display" % m.get("id")
    m.set("matchtype", alignment_type)
    w.update(m)

