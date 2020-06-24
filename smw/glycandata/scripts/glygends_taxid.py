#!/bin/env python27

import sys, urllib, json
from collections import defaultdict
from optparse import OptionParser
import csv

datasets_data = """
GLYDS000142	human_proteoform_glycosylation_sites_harvard.csv	9606	HarvardU
GLYDS000335	hcv1a_proteoform_glycosylation_sites_literature.csv	11108
GLYDS000479	sarscov2_proteoform_glycosylation_sites_unicarbkb.csv	2697049
GLYDS000480	human_proteoform_glycosylation_sites_gptwiki.csv	9606	GPTwiki	saccharide
GLYDS000481	human_proteoform_glycosylation_sites_literature_mining.csv	9606
GLYDS000492	mouse_proteoform_glycosylation_sites_literature_mining.csv	10090
GLYDS000493	rat_proteoform_glycosylation_sites_literature_mining.csv	10116
"""

datasets = dict()
for l in datasets_data.splitlines():
    if not l.strip():
	continue
    sl = l.split('\t')
    sl[2] = int(sl[2])
    datasets[sl[0]] = dict(dsid=sl[0],filename=sl[1],taxid=sl[2])
    if len(sl) >= 4:
	datasets[sl[0]]['source'] = sl[3]
    if len(sl) >= 5:
	datasets[sl[0]]['sourceid'] = sl[4]

seen = set()
for ds in sorted(datasets.values(),key=lambda d: d['dsid']):
  url = "https://data.glygen.org/ln2wwwdata/reviewed/%s"%(ds['filename'],)
  for row in csv.DictReader(urllib.urlopen(url)):
    accs = []
    for acckey in ("saccharide",):
	if row.get(acckey) != None:
	    accs.append(row.get(acckey))
    accs = filter(lambda s: bool(s),accs)
    # if len(accs) == 0:
    # 	print ds['dsid'],row
    if len(accs) == 0:
	continue
    taxid = ds['taxid']
    if 'source' in ds:
        source = ds['source']
    else:
	source = ds['dsid']
    for acc in accs:
        sourceid = None
	if 'sourceid' in ds:
	    sourceid = row.get(ds['sourceid'])
        if (acc,taxid,source,sourceid) in seen:
	    continue
        seen.add((acc,taxid,source,sourceid))
        if sourceid:
	    print "\t".join(map(str,[acc,taxid,source,sourceid]))
        else:
	    print "\t".join(map(str,[acc,taxid,source]))
	break
