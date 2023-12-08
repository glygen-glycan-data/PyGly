#!/bin/env python2

import os, ssl, time

import sys, urllib, json
from collections import defaultdict
from optparse import OptionParser
import csv

import findpygly
from pygly.GlycanResource import GlyGenSourceFile, GlyGenDataset, GlyGenSourceFile

ggds = GlyGenDataset(verbose=True)
for acc,taxid,dsid in ggds.alltaxa():
    print "\t".join(map(str,[acc,taxid,"GlyGen",dsid]))

ggsf = GlyGenSourceFile(verbose=False)
acc2gtc = defaultdict(set)
for row in ggsf.allsourcegtc():
    if row[1] == "-":
        continue
    if row[0] != row[1]:
        acc2gtc[row[0]].add(row[1])

seen = set()
for origid,gtcacc,taxid,source,sourceid in ggsf.allsourcetaxa():
    if source == "GlyConnect":
        if sourceid.startswith("S"):
            source = "GlyConnectStructure"
        else:
            source = "GlyConnectComposition"
        accs = acc2gtc[sourceid]
        sourceid = sourceid[1:]
    elif sourceid in acc2gtc:
        accs = acc2gtc[origid]
    else:
        accs = [gtcacc]
    for acc in accs:
        if (acc,taxid,source,sourceid) in seen:
            continue
        seen.add((acc,taxid,source,sourceid))
        # srcout = source2out.get(source,"GlyGen")
        # if srcout == "GlyGen":
        #     sourceid = sourceid2out[(source,section)]
        if sourceid != None:
            print "\t".join(map(str,[acc,taxid,source,sourceid]))
        else:
            print "\t".join(map(str,[acc,taxid,source]))
