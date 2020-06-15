#!/bin/env python27

import sys
from getwiki import GPTWiki
import numpy as np
from collections import defaultdict

headers = """
PrecursorMz ProductMz Tr_recalibrated
transition_name
CE LibraryIntensity
transition_group_id
decoy PeptideSequence
ProteinName Annotation
FullUnimodPeptideName
MissedCleavages Replicates NrModifications
PrecursorCharge GroupLabel UniprotID
"""

headers = """
peptide_id
transition_group_id
transition_id
transition_group_name
transition_name
PrecursorMz
PrecursorCharge
ProductMz
ProductCharge
Tr_recalibrated
LibraryIntensity
""".split()

w = GPTWiki()

sampids = set(sys.argv[1:])
spec2sampid = dict()

peps = defaultdict(lambda: defaultdict(dict))

i = 0
for tg in w.itertransgroups():
    spec = w.get(tg.get('spectra'))                                                                                         
    if not spec in spec2sampid:
        spec2sampid[spec] = spec.get('sample')
    if not spec2sampid[spec] in sampids:
	continue
    pep = w.get(tg.get('peptide'))
    pepid = pep.get('id')
    if pep.get('nrt') == None:
	continue
    z1 = tg.get('z1')
    ntrans = len(tg.get('transitions',[]))
    if ntrans == 0:
	continue
    print tg.get('id')
    for trid,trint in tg.get('transitions',[]):
	tr = w.get(trid)
	if trid not in peps[(pepid,z1)]:
            peps[(pepid,z1)][trid]['pepname'] = pep.get('name')
            peps[(pepid,z1)][trid]['nrt'] = pep.get('nrt')
            peps[(pepid,z1)][trid]['label'] = tr.get('label')
            peps[(pepid,z1)][trid]['mz2'] = tr.get('mz2')
            peps[(pepid,z1)][trid]['z2'] = tr.get('z2')
            peps[(pepid,z1)][trid]['mz1'] = tr.get('mz1')
            peps[(pepid,z1)][trid]['z1'] = tr.get('z1')
            peps[(pepid,z1)][trid]['int'] = []
        peps[(pepid,z1)][trid]['int'].append(trint)
    i += 1
    if False and i >= 100:
	break

for chpep in peps:
    for trid in peps[chpep]:
	peps[chpep][trid]['int'] = np.median(peps[chpep][trid]['int'])

print "\t".join(headers)
for chpep in peps:
    bpi = None
    for trid in sorted(peps[chpep],key=lambda trid: peps[chpep][trid]['int'], reverse=True):
	tr = peps[chpep][trid]
	if bpi == None:
	    bpi = tr['int']
	relint = int(round(1000.0*tr['int']/bpi,0))
	row = {}
	row['peptide_id'] = chpep[0]
	row['transition_group_id'] = "%s/%d"%chpep
	row['transition_group_name'] = "%s/%d"%(tr['pepname'],chpep[1])
	row['transition_id'] = trid
	row['transition_name'] = "%s/%d"%(tr['label'],tr['z2'])
	row['PrecursorMz'] = tr['mz1']
	row['PrecursorCharge'] = tr['z1']
	row['ProductMz'] = tr['mz2']
	row['ProductCharge'] = tr['z2']
	row['Tr_recalibrated'] = tr['nrt']
	row['LibraryIntensity'] = relint
	print "\t".join(map(str,map(row.get,headers)))
