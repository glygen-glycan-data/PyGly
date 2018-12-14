#!/bin/env python27
import json

from getwiki import GlycoMotifWiki
w = GlycoMotifWiki()

from pygly.GlyTouCan import GlyTouCan
gtc = GlyTouCan()

brokenImage = json.loads(open("brokenimage.json").read())
for m in w.itermotif():
    if m.get("collection") == "GM" and m.get("glytoucan") not in brokenImage:
        m.set("displayhgv", True)
    else:
        m.set("displayhgv", False)
    w.update(m)
