#!/bin/env python3.12

import os, ssl, time

import sys, urllib, json
from collections import defaultdict
from optparse import OptionParser
import csv

import findpygly
from pygly.GlycanResource import GlyGenSourceFile, GlyGenDataset, GlyGenSourceFile

ggds = GlyGenDataset(verbose=True)
for acc,taxid,dsid in ggds.alltaxa():
    print("\t".join(map(str,[acc,taxid,"GlyGen",dsid])))

ggsf = GlyGenSourceFile(verbose=False)
acc2gtc = defaultdict(lambda: defaultdict(set))
for row in ggsf.allsourcegtc():
    if row[1] == "-":
        continue
    if row[0] and row[0] != row[1]:
        acc2gtc[row[2]][row[0]].add(row[1])

seen = set()
for origid,gtcacc,taxid,source,sourceid in ggsf.allsourcetaxa():
    # print(origid,gtcacc,taxid,source,sourceid)
    accs = None
    if source == "GlyConnect":
        if origid.startswith("S"):
            source = "GlyConnectStructure"
        else:
            source = "GlyConnectComposition"
        accs = acc2gtc['GlyConnect'][origid]
        sourceid = origid[1:]
    elif source == "GlyCosmos":
        accs = [ origid ]
        sourceid = origid
    elif sourceid in acc2gtc[source]:
        accs = acc2gtc[source][origid]
    elif gtcacc not in ("-",None):
        accs = [gtcacc]
    else:
        continue
    for acc in accs:
        if (acc,taxid,source,sourceid) in seen:
            continue
        seen.add((acc,taxid,source,sourceid))
        # srcout = source2out.get(source,"GlyGen")
        # if srcout == "GlyGen":
        #     sourceid = sourceid2out[(source,section)]
        if sourceid != None:
            print("\t".join(map(str,[acc,taxid,source,sourceid])))
        else:
            print("\t".join(map(str,[acc,taxid,source])))
