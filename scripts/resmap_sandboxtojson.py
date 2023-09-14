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

sandbox = GlycoTreeSandbox()

if len(accs) == 0:
    accs = sandbox.list()

# Loop through the file paths and read the contents of each file
for acc in accs:
    json_filename = os.path.join(out_dir,acc+".json")
    if not os.path.exists(json_filename):
        continue

    structure_dict = json.loads(open(json_filename).read())
    canon_parent_id = {}
    for res in structure_dict['residues']:
        if res.get('parentid'):
            canon_parent_id[str(res['residueid'])] = str(res.get('parentid'))

    enzymes = dict()
    enzymes['__synonyms__'] = dict()

    jdoc = sandbox.glycan(acc)
    for r in jdoc['residues']:
        can_res_index = r.get('canonical_residue_index')
        if not can_res_index:
            continue
        for enz in r['enzymes']:
            up = enz['uniprot']
            gn = enz['gene_name']
            sp = enz['species']
            if sp == 'Homo sapiens':
                org = 'Human'
            if sp == 'Mus musculus':
                org = 'Mouse'
            if up not in enzymes:
                enzymes[up] = ([],[])
            enzymes[up][0].append(str(can_res_index))
            if canon_parent_id.get(str(can_res_index)):
                enzymes[up][1].append("%s-%s"%(canon_parent_id[str(can_res_index)],can_res_index))
            enzymes['__synonyms__'][gn] = up
            enzymes['__synonyms__']["%s.%s"%(org,gn)] = up
    for up in enzymes:
        if up == '__synonyms__':
            continue
        enzymes[up] = sorted(set(enzymes[up][0]),key=int) + \
                      sorted(set(enzymes[up][1]),key=lambda t: tuple(map(int,t.split('-'))))
    
    structure_dict['annotations']['Enzyme'] = enzymes
    print("%s.json -> %s.json"%(acc,acc))
    with open(json_filename, "w") as json_file:
        json.dump(structure_dict, json_file, indent=4, sort_keys=True)
