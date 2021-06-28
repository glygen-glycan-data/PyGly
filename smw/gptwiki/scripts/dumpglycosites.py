#!/bin/env python2

import sys
from collections import defaultdict
from getwiki import GPTWiki

w=GPTWiki()

tla = {'N': 'Asn'}

sampleids = set(sys.argv[1:])

prsites = defaultdict(dict)

for i,tg in enumerate(w.itertransgroups()):
    pep = w.get(tg.get('peptide'))
    # if not pep.get('nrt'):
    #     continue
    glyid = pep.get('glycan')[0][0]
    gly = w.get(glyid)
    name = ""
    for n in gly.get('name',[]):
	if 'Fuc' in n:
	    continue
	if 'HexNAc' in n:
	    name = n
	    break
    cls = gly.get('class')
    if len(cls) != 1:
	continue
    cls = cls[0]
    spec = w.get(tg.get('spectra'))
    if not spec.get('sample') in sampleids:
        continue
    aligns = pep.get('alignments',[])
    if len(aligns) != 1:
        continue
    for al in aligns:
        pr = al.get('protein')
        site = al.get('prsites')
	if glyid not in prsites[(pr,site[0],int(site[1:]))]:
	    prsites[(pr,site[0],int(site[1:]))][glyid] = dict(aa=tla[site[0]],name=name,cls=cls,pep=set([pep.get('id')]))
	else:
            prsites[(pr,site[0],int(site[1:]))][glyid]['pep'].add(pep.get('id'))

headers = ['UniProt','Site','AminoAcid', 'GlyTouCan', 'Composition', 'GlycoType','Glycopeptides','SiteLink']

print "\t".join(headers)
for (pracc,aa,site),glycans in sorted(prsites.items()):
    for glyid,data in sorted(glycans.items()):
	pepstr = ",".join(sorted(data['pep']))
        row = {'UniProt': pracc, 'Site': site,
               'AminoAcid': data['aa'], 'GlyTouCan': glyid,
	       'Composition': data['name'], 'GlycoType': data['cls'], 
               'Glycopeptides': pepstr, 'SiteLink': "%s@%s%s"%(pracc,aa,site)}
        print "\t".join(map(str,map(row.get,headers)))
