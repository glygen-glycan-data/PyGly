#!/bin/env python27

from getwiki import GlycoMotifWiki, Enzyme
import sys, re, glob, json

w = GlycoMotifWiki()

allenz = dict()
tree = dict()

for jf in glob.glob(sys.argv[1]+'/*.json'):
    enzdata = json.loads(open(jf).read())
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

for jf in glob.glob(sys.argv[1]+'/*.json'):
    enzdata = json.loads(open(jf).read())
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
            allenz[e['gene_name']] = dict(genename=e['gene_name'],
                                          species=('Human' if e['species'] == 'Homo sapiens' else 'Mouse'),
                                          uniprot=e['uniprot'],action=set([action]))

for ed in allenz.values():
    if 'action' in ed:
      d = []
      for i,a in enumerate(sorted(ed['action'])):
        if a[3]:
            d.append('%s %s %s %s to %s'%("Attach" if i == 0 else "attach",a[1],a[2],a[0],a[3]))
        else:
            d.append('%s %s to Asn'%("Attach" if i == 0 else "attach",a[0]))
      ed['description'] = "; ".join(d) + "."
      del ed['action']
    else:
      ed['description'] = ""
    enz = w.get(ed['genename'])
    enz.update(**ed)
    if not enz:
        enz = Enzyme(**ed)
    enz.delete('action')
    if w.put(enz):
        print("Enzyme %s updated."%(ed['genename'],))
