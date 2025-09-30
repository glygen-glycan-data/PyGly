#!/bin/env python3.12

from getwiki import GlycoMotifWiki, GlyTouCanMotif, Class
import sys

w = GlycoMotifWiki()

from clseng import ClassifierEngine

for cl in ClassifierEngine().classifiers():
    theid = cl._id
    wcl = w.get(cl._id)
    if not wcl:
        wcl = Class(id=cl._id)
    wcl.update(hasmotif=cl._motifs,nothasmotif=cl._exceptions,name=cl._class[-1],
               level=("Type" if len(cl._class) == 1 else "Subtype"),
               type=(cl._class[0] if len(cl._class) > 1 else None))
    if w.put(wcl):
        print(wcl.astemplate())
        # sys.exit()

