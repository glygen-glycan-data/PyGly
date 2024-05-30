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

accs = sys.argv[1:]

sandbox = GlycoTreeSandbox(local=True,delaytime=3)

def glycansiter(sandbox,accs):
    if len(accs) == 0:
        for jd in sandbox.allglycans(blocksize=50):
            yield jd
    else:
        for acc in accs:
            yield sandbox.glycan(acc)

headers = """
    glytoucan_ac residue_name residue_id uniprot gene_name gene_id parent_residue_id enzyme_type species
""".split()

print("\t".join(headers))
# Loop through the file paths and read the contents of each file
for sbjdoc in glycansiter(sandbox,accs):
    acc = sbjdoc['accession']
    # print(json.dumps(sbjdoc))
    for res in sbjdoc['residues']:
        if len(res['enzymes']) == 0:
            continue
        for enz in res['enzymes']:
            data = dict(glytoucan_ac=acc,residue_name=res['residue_name'],residue_id=res['residue_id'],
                        uniprot=enz['uniprot'],gene_name=enz['gene_name'],gene_id=enz['gene_id'],
                        parent_residue_id=res['parent_id'],enzyme_type=enz['type'],species=enz['species'])
            print("\t".join(map(lambda h: data[h],headers)))

