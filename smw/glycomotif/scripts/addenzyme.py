#!/bin/env python27

from getwiki import GlycoMotifWiki, Enzyme
import sys, re, glob, json
import findpygly
from pygly.GlycanResource import GlycoTreeSandbox, GlycoTreeSandboxDev

# gts = GlycoTreeSandbox()
gts = GlycoTreeSandboxDev()

allenz = dict()
tree = dict()

for gtsselector in ['mapped_N','mapped_O']:
  for enzdata in gts.allglycans(gtsselector):
    print(enzdata['accession'])
    for r in enzdata["residues"]:
        if r.get('glycotree',"none") == "none":
            continue
        glycotree = r['glycotree']
        residue_id = r["residue_id"]
        name = r["name"]
        anomer = r["anomer"]
        if anomer == "a":
            anomer = "alpha"
        elif anomer == "b":
            anomer = "beta"
        site = r["site"]
        parent_id = r["parent_id"]
        tree[residue_id] = dict(id=residue_id, name=name, anomer=anomer, site=site, 
                                parent_id=parent_id, glycotree=r['glycotree'])

for r in tree.values():
    if r['parent_id'] == '0':
        continue
    p = tree[r['parent_id']]
    ch = p.get('child',set())
    ch.add(r['id'])
    p['child'] = ch

for gtsselector in ['mapped_N','mapped_O']:
  for enzdata in gts.allglycans(gtsselector):
    print(enzdata['accession'])
    for r in enzdata["residues"]:
        if r.get('glycotree',"none") == "none":
            continue
        rid = r['residue_id']
        pid = r['parent_id']
        action = (tree[rid]['name'],tree[rid]['anomer'],tree[rid]['site'],tree[pid]['name'] if pid != "0" else "")
        for e in r['enzymes']:
            if e['gene_name'] in allenz:
                allenz[e['gene_name']]['action'].add(action)
                continue
            allenz[e['gene_name']] = dict(genename=e['gene_name'],glycotree=gtsselector.split('_')[-1],
                                          species=('Human' if e['species'] == 'Homo sapiens' else 'Mouse'),
                                          uniprot=e['uniprot'],action=set([action]))

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
    if 'glycotree' in ed:
        del ed['glycotree']
    enz = w.get(ed['genename'])
    if not enz:
        enz = Enzyme(**ed)
    else:
        enz.update(**ed)
    enz.delete('action')
    enz.delete('tree')
    if w.put(enz):
        print("Enzyme %s updated."%(ed['genename'],))
