#!/bin/env python3.12

import sys
from operator import itemgetter
from collections import defaultdict

from getwiki import GlycoMotifWiki
w = GlycoMotifWiki()

import findpygly

from pygly.GlycanFormatter import WURCS20Format, GlycanParseError
from pygly.GlycanResource import GlyTouCan, GlyTouCanNoCache, GlyTouCanNoPrefetch, GlyCosmosNoCache

debug = False
gtccache = True
gtc = GlyTouCan(usecache=gtccache)

def iterglycan():
    global debug
    if len(sys.argv) > 1:
        debug = True
        for acc in sys.argv[1:]:
            yield acc
    else:
        archived = set()
        gco = GlyCosmosNoCache()
        for acc in gco.archived():
            archived.add(acc["accession"])
        gtc = GlyTouCan(usecache=gtccache)
        for acc,seq,fmt in gtc.allseq(format='wurcs'):
            if acc in archived:
                continue
            yield acc

allacc = set(iterglycan())

classmotif = dict()
classmotifacc = set()
for acc,gtcacc,altype in w.itermotifgtcalign(regex='^GGM.001'):
    assert int(acc.split('.')[1]) >= 1000
    classmotif[(gtcacc,altype)] = acc
    classmotifacc.add(gtcacc)

struct2motif = defaultdict(lambda: defaultdict(set))
headers = "Core_Inclusive  Substructure_Inclusive  Whole_Inclusive Non_Red_Inclusive       Core_Strict     Substructure_Strict     Whole_Strict    Non_Red_Strict".split()
altypemap = dict(zip(["Core","Substructure","Whole-Glycan","Nonreducing-End"],
                     ["Core","Substructure","Whole","Non_Red"]))
for l in sys.stdin:
    if not l.startswith('G'):
        continue
    sl = l.split()
    mgtcacc = sl[0]
    if mgtcacc not in classmotifacc:
        continue
    sgtcacc = sl[1]
    if debug and sgtcacc not in allacc:
        continue
    data = dict(zip(headers,sl[2:]))
    for alt in ("Core","Substructure","Whole-Glycan","Nonreducing-End"):
        if (mgtcacc,alt) not in classmotif:
            continue
        if data[altypemap[alt]+'_Strict'][0] == "Y":
            struct2motif[sgtcacc]['Strict'].add(classmotif[(mgtcacc,alt)])
        if data[altypemap[alt]+'_Inclusive'][0] == "Y":
            struct2motif[sgtcacc]['Loose'].add(classmotif[(mgtcacc,alt)])

from clseng import ClassifierEngine
classifier = ClassifierEngine(alignments=struct2motif,glytoucan=gtc,verbose=debug)

acc2type = defaultdict(set)
for acc in sorted(allacc):

    for clid,asn,motif in classifier.assign(acc):
        acc2type[acc].add((clid,asn,motif))

for acc in sorted(allacc):

    subtypes = defaultdict(set)
    for clid,asn,motif in acc2type[acc]:
        subtypes[asn].add(clid)

    subtypecnt = 0
    stbytcnt = defaultdict(int)
    typecnt = 0
    rows = defaultdict(set)
    for st in sorted(subtypes,key=lambda t: len(t)):
        for sid in sorted(subtypes[st]):
            for i,sti in enumerate(st):
                stt = ("Type" if i == 0 else "Subtype")
                rows[(stt,sti)].add(sid)
    for key,stis in rows.items():
        if key[0] == "Type":
            print("\t".join([acc,key[0],key[1],",".join(sorted(stis))]))
    for key,stis in rows.items():
        if key[0] == "Subtype":
            print("\t".join([acc,key[0],key[1],",".join(sorted(stis))]))
