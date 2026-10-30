#!/bin/env python3.12

import os, ssl, time

import sys, urllib, json
from collections import defaultdict
from optparse import OptionParser
import csv

import findpygly
from pygly.GlycanResource import GlyGenSourceFile, GlyGenDataset, GlyGenSourceFile

ggsf = GlyGenSourceFile(verbose=False)
gtc2acc = defaultdict(lambda: defaultdict(set))
for row in ggsf.allsourcegtc():
    if row[1] == "-":
        continue
    if row[0] and row[0] != row[1]:
        gtc2acc[row[2]][row[1]].add(row[0])

for source in gtc2acc:
    for gtcacc in gtc2acc[source]:
        for acc in gtc2acc[source][gtcacc]:
            if source == "GlyConnect":
                if acc.startswith("S"):
                    print("\t".join(map(str,[acc[1:],gtcacc,"GlyConnectStructure"])))
                else:
                    print("\t".join(map(str,[acc[1:],gtcacc,"GlyConnectComposition"])))
            elif source == "UniCarbKB":
                pass
            else:        
                print("\t".join(map(str,[acc,gtcacc,source])))
