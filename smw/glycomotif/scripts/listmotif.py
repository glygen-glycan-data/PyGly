#!/bin/env python3.12

from getwiki import GlycoMotifWiki, GlyTouCanMotif, Class
import sys

w = GlycoMotifWiki()

from clseng import ClassifierEngine

ce = ClassifierEngine()

for cl in ce.classifiers():
    theid = cl._id
    wcl = w.get(cl._id)
    if not wcl:
        wcl = Class(id=cl._id)
    wcl.update(hasmotif=cl._motifs,nothasmotif=cl._exceptions,name=cl._class[-1],
               level=("Type" if len(cl._class) == 1 else "Subtype"))
    if w.put(wcl):
        print(cl._id)

sys.exit()
for cl in w.iterpages(regex=r'^GGC\.'):
    print(w.get(cl.name))
sys.exit()

for r in w.itermotifgtcalign(regex=r'^GGM\.001'):
    print(r)
sys.exit()

for m in w.itermotif(regex='GGM.001'):
    keywords = set(m.get('keyword',[]))
    keywords.add("Glycan classification")
    m.set('keyword',keywords)
    if w.put(m):
        print(m.get('id'))
