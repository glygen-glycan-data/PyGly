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
from pygly.alignment import GlycanEqual, GlycanTopoEqual, GlycanImageEqual, GlycanCompImageEqual
from pygly.GlycanFormatter import *
from pygly.MonoFormatter import *
from pygly.GlycanBuilderSVGParser import GlycanBuilderSVG

ctable = ResidueCompositionTable()
iupacSym = IUPACSym()

ip = IUPACLinearFormat()
wurcs_dir = sys.argv[1]
svg_dir = sys.argv[2]

accs = sys.argv[3:]

wp = WURCS20Format()
gctp = GlycoCTFormat()
jcr = JSONCanonResidue()
sp = GlycanBuilderSVG()
glyeq = GlycanEqual()
glyimeq = GlycanImageEqual()
glycompimeq = GlycanCompImageEqual()

# Loop through the file paths and read the contents of each file
for path in sorted(os.listdir(wurcs_dir)):

    if not path.endswith('.txt'):
        continue

    acc = os.path.splitext(os.path.basename(path))[0]
    if len(accs) != 0 and acc not in accs:
        continue
    print(acc)

    # We must have an SVG file...
    svg_filename = os.path.join(svg_dir,acc+".svg")
    if not os.path.exists(svg_filename):
        continue
    
    canon_seq = open(os.path.join(wurcs_dir,path), "r").read().strip()
    canon_seqhash = hashlib.md5(canon_seq).hexdigest().lower()

    try: 
        canon_gly = wp.toGlycan(canon_seq)
    except GlycanParseError:
        # traceback.print_exc(file=sys.stdout)
        continue

    svg_seq = open(svg_filename, "r").read().strip()
    svg_seqhash = hashlib.md5(svg_seq).hexdigest().lower()

    try:
        svg_gly = sp.toGlycan(svg_seq)
    except GlycanParseError:
        # traceback.print_exc(file=sys.stdout)
        continue

    canonres_data = {}
    iupac_annotations = defaultdict(list)
    for m in canon_gly.all_nodes(undet_subst=True):
        data = json.loads(jcr.toStr(m))
        canonres_data[str(m.id())] = data

    for mid,iupacsym,isaggr in canon_gly.iupac_items(canon_gly.all_nodes(undet_subst=True)):
        iupac_annotations[iupacsym].append(str(mid))

    canonres_data = sorted(canonres_data.values(),key=lambda item: int(item['residueid']))

    svg_idmap = []
    if canon_gly.has_root() and svg_gly.has_root():
        if not glyimeq.eq(svg_gly,canon_gly,idmap=svg_idmap):
            continue
    elif not canon_gly.has_root() and not svg_gly.has_root():
        if not glycompimeq.eq(svg_gly,canon_gly,idmap=svg_idmap):
            continue
    else:
        continue

    svg_idmapids = [ (t[0].external_descriptor_id(),t[1].id()) for t in svg_idmap ]

    svg_idmap_dict = defaultdict(list)
    for svg_id,canon_id in svg_idmapids:
        svg_idmap_dict[str(canon_id)].append(svg_id)
        
    for l in canon_gly.all_links():
        for parent_svgid in svg_idmap_dict[str(l.parent().id())]:
            for child_svgid in svg_idmap_dict[str(l.child().id())]:
                svgidbase,parent_svgid1 = parent_svgid.rsplit(':',1)
                svgidbase = svgidbase.split('-',1)[1]
                child_svgid1 = child_svgid.rsplit(':',1)[1]
                link_id = str(l.parent().id()) + "-" + str(l.child().id())
                svg_link_id = "l-1:" + str(parent_svgid1) + ","+ str(child_svgid1)
                svg_idmap_dict[link_id].append(svg_link_id)

    structure_dict = {}
    structure_dict['canonical_sequence_accession'] = acc
    structure_dict['canonical_sequence_md5'] = canon_seqhash
    
    structure_dict['residues'] = canonres_data
    
    structure_dict['residuemaps'] = {}
    structure_dict['residuemaps'][svg_seqhash+':SVG'] = svg_idmap_dict
    
    structure_dict['annotations'] = {}
    structure_dict['annotations']['IUPAC'] = iupac_annotations
    
    jsonfilename = acc + ".json"
    with open(jsonfilename, "w") as json_file:
        json.dump(structure_dict, json_file, indent=4, sort_keys=True)
