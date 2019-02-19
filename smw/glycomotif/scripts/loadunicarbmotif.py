#!/bin/env python27

import sys, traceback
from getwiki import GlycoMotifWiki, UniCarbMotif
w = GlycoMotifWiki()

from pygly.GlyTouCan import GlyTouCan
gtc = GlyTouCan()

from gtccache import GlyTouCanCache
gtccache = GlyTouCanCache()

from pygly.GlycanFormatter import GlycoCTFormat, IUPACParserExtended1
gparser = GlycoCTFormat()
imparser = IUPACParserExtended1()

from dataset import XLSXFileTable
rows = XLSXFileTable(sys.argv[1])

possibleaglycon = ["Cer", "R", "Ser/Thr"]
reaglycon = ["Ser/Thr", "Cer", "Other"]
current = set()
for r in rows:
    id = r["ID"]
    name = r["Name"]
    iupacseq = r["IUPAC"]

    accession = "%06d" % id
    if not iupacseq:
        continue

    if "n-" in iupacseq:
        print "%s Contains repeats" % accession
        continue

    try:
        gobj = imparser.toGlycan(iupacseq)
        glycoct = gparser.toStr(gobj)
        wurcs = gtc.glycoct2wurcs(glycoct)
        glytoucan = gtc.register(wurcs)[0]

    except:
        print "%s is not able to load to Glycomotif" % accession
        print iupacseq
        continue

    aglycon = None
    redend = None
    for agly in possibleaglycon:
        if agly in iupacseq:
            aglycon = agly
            if aglycon in reaglycon:
                redend = True
            else:
                redend = False
            break

    if False:
        print accession, name
        print iupacseq
        print glycoct
        print glytoucan
        print aglycon
        print redend
        print "\n\n\n\n"

    motif = UniCarbMotif(accession=accession, name=name, glytoucan=glytoucan, redend=redend, aglycon=aglycon,
                    glycoct=glycoct, wurcs=wurcs)
    if not motif.get('glycoct'):
        motif.set('glycoct', glycoct)
    if w.put(motif):
        print >> sys.stderr, accession
    current.add(accession)

for m in w.itermotif():
    if UniCarbMotif.id not in m.get("pagename"):
        continue
    if m.get('accession') not in current:
        print >> sys.stderr, "Deleting:", m.get('pagename')
        w.delete(m.get('pagename'))