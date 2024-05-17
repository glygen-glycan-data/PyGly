#!/bin/env python27

from getwiki import GlycoMotifWiki, Enzyme
import sys, re, glob, json
import findpygly
from pygly.GlycanResource import GlycoTreeSandbox, GlycoTreeSandboxDev

gts = GlycoTreeSandbox()
# gts = GlycoTreeSandboxDev()

form_name_dict = dict(map(lambda s: tuple(s.split()),filter(lambda s: s.strip(),"""
Fucp Fuc
Fucx Fuc
Galf Gal
Galp Gal
GalpNAc GalNAc
GalxNAc GalNAc
GlcAp GlcA
Glcp Glc
Glcx Glc
GlcpNAc GlcNAc
GlcxNAc GlcNAc
KDN KDN
Manp Man
Manx Man
ManpNAc ManNAc
NeupNAc NeuAc
NeupNGc NeuGc
phosphate phosphate
sulfate sulfate
Xylf Xyl
Xylp Xyl
""".splitlines())))

allenz = dict()
tree = dict()

for r in gts.enzymes():
    glycotree = r['residue_id'][0]
    residue_id = r["residue_id"]
    name = form_name_dict[r["form_name"]]
    anomer = r["anomer"]
    if anomer == "a":
        anomer = "alpha"
    elif anomer == "b":
        anomer = "beta"
    site = r["site"]
    parent_id = r["parent_id"]
    if parent_id == 'no_id':
        parent_id = None
    tree[residue_id] = dict(id=residue_id, name=name, anomer=anomer, site=site, 
                            parent_id=parent_id, glycotree=glycotree)
    if r.get('gene_name'):
        if r['gene_name'] not in allenz:
            allenz[r['gene_name']] = dict(genename=r['gene_name'],glycotree=glycotree,
                                          species=('Human' if r['species'] == 'Homo sapiens' else 'Mouse'),
                                          uniprot=r['uniprot'],residues=set([residue_id]))
        else:
            allenz[r['gene_name']]['residues'].add(residue_id)

for rid in tree:
    pid = tree[rid]['parent_id']
    action = (tree[rid]['name'],tree[rid]['anomer'],tree[rid]['site'],tree[pid]['name'] if pid else None)
    tree[rid]['action'] = action

for gn in allenz:
    allenz[gn]['action'] = set()
    for rid in allenz[gn]['residues']:
        allenz[gn]['action'].add(tree[rid]['action'])
    del allenz[gn]['residues']

w = GlycoMotifWiki()

for ed in allenz.values():
    if 'action' in ed:
      d = []
      for i,a in enumerate(sorted(ed['action'])):
        if a[3]:
            if a[0] not in ('sulfate',):
                d.append('%s %s %s %s to %s'%("Attach" if i == 0 else "attach",a[1],a[2],a[0],a[3]))
            else:
                d.append('%s %s linked %s to %s'%("Attach" if i == 0 else "attach",a[2],a[0],a[3]))
        else:
            if ed['glycotree'] == 'N':
                d.append('%s %s to Asn'%("Attach" if i == 0 else "attach",a[0]))
            elif ed['glycotree'] == 'O':
                d.append('%s %s to Ser/Thr'%("Attach" if i == 0 else "attach",a[0]))
      ed['description'] = "; ".join(d) + "."
      del ed['action']
    else:
      ed['description'] = ""
    if 'action' in ed:
        del ed['action']
    if 'residues' in ed:
        del ed['residues']
    if 'glycotree' in ed:
        del ed['glycotree']
    enz = w.get(ed['genename'])
    if not enz:
        enz = Enzyme(**ed)
    else:
        enz.update(**ed)
    enz.delete('action')
    enz.delete('tree')
    enz.delete('residues')
    if w.put(enz):
        print("Enzyme %s updated."%(ed['genename'],))
