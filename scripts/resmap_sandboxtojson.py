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

sandbox = GlycoTreeSandbox(local=True,delaytime=3)

def glycansiter(sandbox,accs):
    if len(accs) == 0:
        for jd in sandbox.allglycans(blocksize=50):
            yield jd
    else:
        for acc in accs:
            yield sandbox.glycan(acc)

# Loop through the file paths and read the contents of each file
for sbjdoc in glycansiter(sandbox,accs):
    acc = sbjdoc['accession']

    json_filename = os.path.join(out_dir,acc+".json")
    if not os.path.exists(json_filename):
        continue

    try:
        structure_dict = json.loads(open(json_filename).read())
    except ValueError:
        raise RuntimeError("Bad JSON format: "+json_filename)
    canon_parent_id = {}
    for res in structure_dict['residues']:
        if res.get('parentid'):
            canon_parent_id[str(res['residueid'])] = str(res.get('parentid'))

    enzymes = dict()
    enzymes['__synonyms__'] = dict()

    enzymesnc = dict()
    enzymesnc['__synonyms__'] = dict()

    for r in sbjdoc['residues']:
        can_res_index = r.get('canonical_residue_index')
        if not can_res_index:
            continue
        for enz in r['enzymes']:
            up = enz['uniprot']
            gn = enz['gene_name']
            sp = enz['species']
            rv = enz.get('rule_violations',[])
            if sp == 'Homo sapiens':
                org = 'Human'
            if sp == 'Mus musculus':
                org = 'Mouse'
            if up not in enzymes:
                enzymes[up] = ([],[])
            enzymes[up][0].append(can_res_index)
            if canon_parent_id.get(can_res_index):
                enzymes[up][1].append("%s-%s"%(canon_parent_id[can_res_index],can_res_index))
            enzymes['__synonyms__'][gn] = up
            enzymes['__synonyms__']["%s.%s"%(org,gn)] = up
            if len(rv) == 0:
                if up not in enzymesnc:
                    enzymesnc[up] = ([],[])
                enzymesnc[up][0].append(can_res_index)
                if canon_parent_id.get(can_res_index):
                    enzymesnc[up][1].append("%s-%s"%(canon_parent_id[can_res_index],can_res_index))
                enzymesnc['__synonyms__'][gn] = up
                enzymesnc['__synonyms__']["%s.%s"%(org,gn)] = up

    for up in enzymes:
        if up == '__synonyms__':
            continue
        enzymes[up] = sorted(set(enzymes[up][0]),key=float) + \
                      sorted(set(enzymes[up][1]),key=lambda t: tuple(map(float,t.split('-'))))

    for up in enzymesnc:
        if up == '__synonyms__':
            continue
        enzymesnc[up] = sorted(set(enzymesnc[up][0]),key=float) + \
                        sorted(set(enzymesnc[up][1]),key=lambda t: tuple(map(float,t.split('-'))))

    structure_dict['annotations']['Enzyme'] = enzymes
    structure_dict['annotations']['EnzymeNoCaveat'] = enzymesnc
    print("%s.json -> %s.json"%(acc,acc))
    with open(json_filename, "w") as json_file:
        json.dump(structure_dict, json_file, indent=4, sort_keys=True)
