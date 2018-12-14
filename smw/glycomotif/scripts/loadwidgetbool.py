#!/bin/env python27
import json

from getwiki import GlycoMotifWiki
w = GlycoMotifWiki()

brokenImage = json.loads(open("brokenimage.json").read())
for m in w.itermotif():
    if m.get("collection") == "GM" and m.get("glytoucan") not in brokenImage:
        print m.get("id")
        m.set("displayhgv", True)
        w.update(m)
    else:
        # continue
        m.set("displayhgv", False)
        w.update(m)