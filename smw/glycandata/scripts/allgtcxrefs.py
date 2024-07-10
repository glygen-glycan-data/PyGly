#!/bin/env python3

import sys, time, traceback
from collections import defaultdict

import findpygly
from pygly.GlycanResource import GlyTouCan
from pygly.GlycanResource import GlyCosmos

gtc = GlyTouCan(verbose=False,usecache=False)
gco = GlyCosmos(verbose=False,usecache=False)

allgco = set(gco.allaccessions())

archived = set(map(lambda d: d['accession'],gco.archived()))
print("%d accessions archived."%(len(archived),))

print("\t".join(["GlyTouCanAccession","Resource","EntryID"]));

def genrefs():
    for gtcacc,resource,entry in gtc.allcrossrefs():
        yield [gtcacc,resource,entry]
    for gtcacc in allgco:
        yield [gtcacc,"glycosmos",gtcacc]

for gtcacc,resource,entry in sorted(genrefs()):
    if gtcacc in archived:
        continue
    if resource in ('glygen',):
        continue
    print("\t".join([gtcacc,resource,entry]))
