#!/bin/env python3

import json
import glob 
import io
import os
import requests
import sys, glob, hashlib, os, os.path, traceback
import findpygly
from pygly.Glycan import *
from pygly.GlycanResource import *

accs = sys.argv[1:]

sandbox = GlycoTreeSandboxDev(local=True,delaytime=3)

def glycansiter(sandbox,accs):
    if len(accs) == 0:
        for jd in sandbox.allglycans(blocksize=50):
            yield jd
    else:
        for acc in accs:
            yield sandbox.glycan(acc)

headers = """
    glytoucan_ac residue_name residue_id canonical_residue_index uniprot gene_name gene_id parent_residue_id enzyme_type species caveat
""".split()

print("\t".join(headers))
# Loop through the file paths and read the contents of each file
for sbjdoc in glycansiter(sandbox,accs):
    # print(json.dumps(sbjdoc))
    if 'rule_violations' not in sbjdoc:
        continue
    acc = sbjdoc['accession']
    residues = dict()
    enzymes = defaultdict(list)
    for res in sbjdoc['residues']:
        residues[res['residue_id']] = copy.copy(res)
        for enz in res.get('enzymes',[]):
            enzymes[res['residue_id'],enz['uniprot']] = enz
    for rule in sbjdoc['rule_violations']:
        res = residues[rule['focus']]
        if rule['enzyme']:
            enz = enzymes[rule['focus'],rule['enzyme']]
            data = dict(glytoucan_ac=acc,residue_name=res['residue_name'],residue_id=res['residue_id'],
                        canonical_residue_index=res['canonical_residue_index'],
                        uniprot=enz['uniprot'],gene_name=enz['gene_name'],gene_id=enz['gene_id'],
                        parent_residue_id=res['parent_id'],enzyme_type=enz['type'],species=enz['species'])
        else:
            data = dict(glytoucan_ac=acc,residue_name=res['residue_name'],residue_id=res['residue_id'],
                        canonical_residue_index=res['canonical_residue_index'],
                        parent_residue_id=res['parent_id'])
        data['caveat'] = rule['assertion'].replace("&alpha;","a").replace("&beta;","b")
        print("\t".join(map(lambda h: data.get(h,""),headers)))
