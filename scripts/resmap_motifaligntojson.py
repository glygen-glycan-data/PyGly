#!/bin/env python2

import json
import glob 
import io
import os
import requests
import sys, glob, hashlib, os, os.path, traceback
import findpygly
from pygly.Glycan import *
from pygly.GlycanResource import *

if os.path.isdir(sys.argv[1]):
    out_dir = sys.argv[1]
    sys.argv.pop(1)
else:
    out_dir = "."

accs = sys.argv[1:]

gm = GlycoMotifNoCache(verbose=False,local=True)
colls = list(gm.collections())

def getmotifs(acc):
    if len(accs) > 0:
        for r in gm.getmotifbystruct(acc):
           yield r
    else:
        for col in colls:
           for r in gm.getmotif(col,acc):
              yield r

def intstrkey(s):
    try:
        return int(s),""
    except ValueError:
        return 1e+20,s

for path in sorted(os.listdir(out_dir)):

    if not path.endswith('.json'):
        continue

    acc = os.path.splitext(path)[0]
    if len(accs) > 0 and acc not in accs:
        continue

    motifs = dict()
    # synonyms = dict()

    for motifacc,altype,isstrict,nodeids,linkids in getmotifs(acc):
        # print(1,motifacc, altype, isstrict, nodeids, linkids)
        if motifacc not in motifs:
            motifs[motifacc] = [set(filter(None,nodeids.split(','))),
                                set(filter(None,linkids.split(',')))]
        else:
            motifs[motifacc][0].update(filter(None,nodeids.split(',')))
            motifs[motifacc][1].update(filter(None,linkids.split(',')))
        # print(2,motifs[motifacc])
    for motifacc in list(motifs):
        # print(3,motifacc,motifs[motifacc])
        # print(4,sorted(motifs[motifacc][0],key=int))
        # print(5,sorted(motifs[motifacc][1],key=lambda t: tuple(map(int,t.split('-')))))
        motifs[motifacc] = sorted(motifs[motifacc][0],key=float) + \
                           sorted(motifs[motifacc][1],key=lambda t: tuple(map(intstrkey,t.split('-'))))
        # print(6,motifs[motifacc])
    
    json_filename = os.path.join(out_dir,path)
    structure_dict = json.loads(open(json_filename).read())
    # motifs['__synonyms__'] = synonyms

    if len(motifs) == 0:
        if 'MotifAlignments' in structure_dict['annotations']:
            del structure_dict['annotations']['MotifAlignments']
    else:
        structure_dict['annotations']['MotifAlignments'] = motifs
    print("%s.json -> %s.json"%(acc,acc))
    with open(json_filename, "w") as json_file:
        json.dump(structure_dict, json_file, indent=4, sort_keys=True)
