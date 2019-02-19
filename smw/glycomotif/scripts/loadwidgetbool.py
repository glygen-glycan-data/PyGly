#!/bin/env python27
from getwiki import GlycoMotifWiki
w = GlycoMotifWiki()

import os
fpath = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../data/brokenimage.txt")

brokenImage = open(fpath).read().split()

for m in w.itermotif():
    if m.get("collection") == "GM" and m.get("glytoucan") not in brokenImage:
        print "%s set to display" % m.get("id")
        m.set("displayhgv", True)
        w.update(m)
    else:
        print "%s set to NOT display" % m.get("id")
        m.set("displayhgv", False)
        w.update(m)
