#!/bin/env python3.12

import sys
from collections import defaultdict

import findpygly
from pygly.GlycanResource import GlyConnect

glyconn = GlyConnect()

species_data = """
9606	Homo sapiens
10090	Mus musculus
10116	Rattus norvegicus
"""

species = dict()
for l in species_data.splitlines():
    sl = map(str.strip,l.split('\t'))
    if len(sl) != 2:
	continue
    species[int(sl[0])] = sl[1]

seen = set()
for sptaxid,spname in species.items():
    for row in glyconn.rowsbyspecies(spname):
	taxid = row.get('taxid')
	if not taxid or sptaxid != int(taxid):
	    continue
	gtcacc = row.get('gtcacc')
	sourceid = row.get('accession')
	if not gtcacc:
	    gtcacc = row.get('compgtcacc')
	    sourceid = row.get('compaccession')
	if not gtcacc or not sourceid:
	    continue
	if (gtcacc,taxid,sourceid) in seen:
	    continue
	seen.add((gtcacc,taxid,sourceid))
	print("\t".join(map(str,[gtcacc,taxid,"GlyConnect",sourceid])))
