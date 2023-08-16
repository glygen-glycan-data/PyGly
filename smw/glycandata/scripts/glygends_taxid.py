#!/bin/env python2

import os, ssl, time

import sys, urllib, json
from collections import defaultdict
from optparse import OptionParser
import csv

import findpygly
from pygly.GlycanResource import GlyGenSourceFile, GlyGenDataset, GlyGenFile

ggds = GlyGenDataset(verbose=True)
for acc,taxid,dsid in ggds.alltaxa():
    print "\t".join(map(str,[acc,taxid,"GlyGen",dsid]))

source2out = dict(gptwiki="GPTwiki",unicarbkb="UniCarbKB",
                  GlyConnectStructure="GlyConnectStructure",
                  GlyConnectComposition="GlyConnectComposition")

ggsf = GlyGenSourceFile(verbose=True)
acc2gtc = defaultdict(set)
for acc,gtc in ggsf.allgtc():
    if acc != gtc:
        acc2gtc[acc].add(gtc)
seen = set()
for sourceid,taxid,source in ggsf.alltaxa():
    if source == "mcw_oglcnac":
        continue
    if source == "glyconnect":
        if sourceid.startswith("S"):
            source = "GlyConnectStructure"
        else:
            source = "GlyConnectComposition"
        accs = acc2gtc[sourceid]
        sourceid = sourceid[1:]
    elif sourceid in acc2gtc:
        accs = acc2gtc[sourceid]
    else:
        accs = [sourceid]
    for acc in accs:
        if (acc,taxid,source,sourceid) in seen:
            continue
        seen.add((acc,taxid,source,sourceid))
        print "\t".join(map(str,[acc,taxid,source2out[source],sourceid]))
