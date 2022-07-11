#!/bin/env python2

import os, ssl, time

if (not os.environ.get('PYTHONHTTPSVERIFY', '') and getattr(ssl, '_create_unverified_context', None)):
    ssl._create_default_https_context = ssl._create_unverified_context

import sys, urllib, json
from collections import defaultdict
from optparse import OptionParser
import csv

datasets_data = """
GLY_000142	human_proteoform_glycosylation_sites_harvard.csv	9606
GLY_000335	hcv1a_proteoform_glycosylation_sites_literature.csv	11108
"""

# GLY_000481	human_proteoform_glycosylation_sites_literature_mining.csv	9606
# GLYDS000492	mouse_proteoform_glycosylation_sites_literature_mining.csv	10090
# GLYDS000493	rat_proteoform_glycosylation_sites_literature_mining.csv	10116
# GLYDS000479	sarscov2_proteoform_glycosylation_sites_unicarbkb.csv	2697049
# GLYDS000480	human_proteoform_glycosylation_sites_gptwiki.csv	9606	GPTwiki	saccharide

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
  url = "https://data.glygen.org/ln2data/releases/data/current/reviewed/%s"%(ds['filename'],)
  # print >>sys.stderr, url
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
	source = 'GlyGen'
    for acc in accs:
        sourceid = None
	if 'sourceid' in ds:
	    sourceid = row.get(ds['sourceid'])
	if sourceid == None and source == 'GlyGen':
	    sourceid = ds['dsid']
        if (acc,taxid,source,sourceid) in seen:
	    continue
        seen.add((acc,taxid,source,sourceid))
        if sourceid:
	    print "\t".join(map(str,[acc,taxid,source,sourceid]))
        else:
	    print "\t".join(map(str,[acc,taxid,source]))
	break
  time.sleep(60)

import findpygly
from pygly.GlycanResource import GlyGenSourceFile

ggsf = GlyGenSourceFile()

acc2gtc = defaultdict(set)
for acc,gtc in ggsf.glyconnect_allgtc():
    acc2gtc[acc].add(gtc)
for sourceid,taxid in ggsf.glyconnect_alltaxa():
    for acc in acc2gtc[sourceid]:
        if sourceid.startswith("S"):
	   source = "GlyConnectStructure"
	else:
	   source = "GlyConnectComposition"
	sourceid = sourceid[1:]
        if (acc,taxid,source,sourceid) in seen:
            continue
        seen.add((acc,taxid,source,sourceid))
        print "\t".join(map(str,[acc,taxid,source,sourceid]))

source = "GPTwiki"
acc2gtc = defaultdict(set)
for acc,gtc in ggsf.gptwiki_allgtc():
    acc2gtc[acc].add(gtc)
for sourceid,taxid in ggsf.gptwiki_alltaxa():
    for acc in acc2gtc[sourceid]:
        if (acc,taxid,source,sourceid) in seen:
            continue
        seen.add((acc,taxid,source,sourceid))
        print "\t".join(map(str,[acc,taxid,source,sourceid]))

