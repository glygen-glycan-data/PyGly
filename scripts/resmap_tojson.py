#!/bin/env python2

from __future__ import print_function

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

verbose = False
if sys.argv[1] == "-v":
    verbose = True
    sys.argv.pop(1)

wurcs_dir = sys.argv[1]
svg_dir = sys.argv[2]
if os.path.isdir(sys.argv[3]):
    out_dir = sys.argv[3]
    sys.argv.pop(3)
else:
    out_dir = "."

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
    filenamebase = acc

    if ":" in acc:
        taskhash,acc = acc.split(':',1)
    elif re.search(r'^G\d{5}[A-Z]{2}$',acc):
        acc = acc
    else:
        acc = None


    # if os.path.exists(os.path.join(out_dir,acc+".json")):
    #     continue

    # We must have an SVG file...
    svg_filename = os.path.join(svg_dir,filenamebase+".svg")
    if not os.path.exists(svg_filename):
        continue
    
    canon_seq = open(os.path.join(wurcs_dir,path)).read().strip()
    canon_seqhash = hashlib.md5(canon_seq.encode('ascii')).hexdigest().lower()

    try: 
        canon_gly = wp.toGlycan(canon_seq)
    except GlycanParseError:
        if verbose:
            print("Error accession:",acc,file=sys.stderr)
            traceback.print_exc(file=sys.stderr)
        continue

    svg_seq = open(svg_filename, "rb").read().strip()
    svg_seqhash = hashlib.md5(svg_seq).hexdigest().lower()

    try:
        svg_gly = sp.toGlycan(svg_seq)
    except GlycanParseError:
        if verbose:
            print("Error accession:",acc,file=sys.stderr)
            traceback.print_exc(file=sys.stderr)
        continue

    svg_seq = svg_seq.decode('utf8')
    svglinkids = set()
    for m in re.findall(r' ID="(l-1:\d+,\d+)" ',svg_seq):
        svglinkids.add(m)
    # somtimes the link is left out but the linkinfo is there...
    for m in re.findall(r' ID="(li-1:\d+,\d+)" ',svg_seq):
        svglinkids.add(m.replace("li-1:","l-1:"))

    canonres_data = {}
    iupac_annotations = defaultdict(list)
    for m in canon_gly.all_nodes(undet_subst=True):
        data = json.loads(jcr.toStr(m))
        canonres_data[m.external_descriptor_id()] = data

    for mid,iupacsym,isaggr in canon_gly.iupac_items(canon_gly.all_nodes(undet_subst=True)):
        iupac_annotations[iupacsym].append(mid[0])
    for iupacsym in iupac_annotations:
        iupac_annotations[iupacsym] = sorted(iupac_annotations[iupacsym],key=float)
    # iupac_synonyms = {}
    # for iupacsym in iupac_annotations:
    #     if iupacsym.lower() != iupacsym:
    #         iupac_synonyms[iupacsym.lower()] = iupacsym
    # iupac_annotations['__synonyms__'] = iupac_synonyms

    canonres_data = sorted(canonres_data.values(),key=lambda item: float(item['residueid']))

    svg_idmap = []
    if canon_gly.has_root() and svg_gly.has_root():
        if not glyimeq.eq(svg_gly,canon_gly,idmap=svg_idmap):
            if verbose:
                print("Error accession:",acc,file=sys.stderr)
                print("Cannot align SVG-based and WURCS-based glycans",file=sys.stderr)
                print(svg_gly.glycoct())
                print(canon_gly.glycoct())
            continue
        svg_idmapids = [ ti for t in svg_idmap for ti in glyimeq.monoidmap(*t)]
    elif not canon_gly.has_root() and not svg_gly.has_root():
        if not glycompimeq.eq(svg_gly,canon_gly,idmap=svg_idmap):
            if verbose:
                print("Error accession:",acc,file=sys.stderr)
                print("Cannot align SVG-based and WURCS-based composition glycans",file=sys.stderr)
            continue
        svg_idmapids = [ ti for t in svg_idmap for ti in glycompimeq.monoidmap(*t)]

    else:
        if verbose:
            print("Error accession:",acc,file=sys.stderr)
            print("Cannot align SVG-based and WURCS-based glycans - composition mismatch",file=sys.stderr)
        continue

    svg_idmap_dict = defaultdict(list)
    for svg_id,canon_id in svg_idmapids:
        svg_idmap_dict[str(canon_id)].extend(svg_id.split(';'))

    for l in canon_gly.all_links(uninstantiated=True):
        parent_svgid =  svg_idmap_dict[str(l.parent().id())][0]
        child_svgid =  svg_idmap_dict[str(l.child().id())][0]
        svgidbase,parent_svgid1 = parent_svgid.rsplit(':',1)
        svgidbase = svgidbase.split('-',1)[1]
        child_svgid1 = child_svgid.rsplit(':',1)[1]
        if l.instantiated():
            link_id = str(l.parent().id()) + "-" + str(l.child().id())
            svg_link_id = "l-1:" + str(parent_svgid1) + ","+ str(child_svgid1)
        else:
            link_id = "-" + str(l.child().id())
            try:
                svg_link_id = [ sli for sli in svglinkids if sli.endswith(","+str(child_svgid1)) ][0]
            except IndexError:
                svg_link_id = None
        if svg_link_id:
            svg_idmap_dict[link_id].append(svg_link_id)

    for canid in svg_idmap_dict:
        svg_idmap_dict[canid] = sorted(set(svg_idmap_dict[canid]))

    if not os.path.exists(os.path.join(out_dir,filenamebase+".json")):
        structure_dict = {}
    else:
        try:
            structure_dict = json.loads(open(os.path.join(out_dir,filenamebase+".json")).read())
        except ValueError:
            raise RuntimeError("Bad JSON format: "+os.path.join(out_dir,filenamebase+".json"))

    if acc:
        structure_dict['canonical_sequence_accession'] = acc
    structure_dict['canonical_sequence'] = canon_seq
    structure_dict['canonical_sequence_md5'] = canon_seqhash
    
    structure_dict['residues'] = canonres_data
    
    structure_dict['residuemap'] = svg_idmap_dict
    structure_dict['svg_sequence_md5'] = svg_seqhash
    
    if 'annotations' not in structure_dict:
        structure_dict['annotations'] = {}
    structure_dict['annotations']['IUPAC'] = iupac_annotations
    
    print("%s.txt,%s.svg -> %s.json"%(filenamebase,filenamebase,filenamebase))
    jsonfilename = os.path.join(out_dir,filenamebase + ".json")
    with open(jsonfilename, "w") as json_file:
        json.dump(structure_dict, json_file, indent=4, sort_keys=True)
