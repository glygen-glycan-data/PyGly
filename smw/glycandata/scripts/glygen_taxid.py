#!/bin/env python3.12

import sys, urllib, json
from collections import defaultdict
from optparse import OptionParser
import csv

parser = OptionParser()
parser.add_option("--glyconnect",type='string',dest="glyconnect",default=None)
parser.add_option("--unicarbkb",type='string',dest="unicarbkb",default=None)

opt,args = parser.parse_args()

def readmapping(filename):
    map = dict()
    for r in csv.reader(open(filename),dialect='excel-tab'):
        if r[0] == "GlyTouCanAccession":
            continue
        gtc = r[0]
        acc = r[1]
        map[r[0]] = r[1]
    return map

gtcmapping = defaultdict(dict)
if opt.glyconnect:
    gtcmapping['GlyConnect'] = readmapping(opt.glyconnect)
if opt.unicarbkb:
    gtcmapping['UniCarbKB'] = readmapping(opt.unicarbkb)
gtcmapping['GPTwiki'] = None

src2src_data = """
UniCarbKB	UniCarbKB
UniCarbKB|Literature	UniCarbKB
Harvard U	HarvardU
GlyConnect	GlyConnect
GPTwiki	GPTwiki
"""

src2src = dict()
for line in src2src_data.splitlines():
    if line.strip() == "":
	continue
    key,value = map(str.strip,line.split('\t',1))
    src2src[key] = value

glygen_glycosylation_data_urls = """
https://raw.githubusercontent.com/GW-HIVE/data_share/master/current/glycosylation_data.csv?token=ABCMDN3JANE7IDOZ2YPXY3K6WLUHA
https://raw.githubusercontent.com/GW-HIVE/data_share/master/current/glycosylation_data_glyconnect.csv?token=ABCMDN3JGW3RFIF5TSHW5XS6WLYXE
""".strip()

import csv

seen = set()
for url in glygen_glycosylation_data_urls.splitlines():
  for row in csv.DictReader(urllib.urlopen(url.strip())):
    accs = []
    for acckey in ("saccharide","glytoucan_ac_structure","glytoucan_ac_composition"):
	if row.get(acckey):
	    accs.append(row.get(acckey))
    taxid = int(row['tax_id'])
    source = src2src[row['source']]
    for acc in accs:
        sourceid = None
        if source in gtcmapping:
	    if gtcmapping[source] == None:
	        sourceid = acc
	    else:
	        sourceid = gtcmapping[source].get(acc,None)
        if (acc,taxid,source,sourceid) in seen:
	    continue
        if sourceid:
	    print("\t".join(map(str,[acc,taxid,source,sourceid])))
        else:
	    print("\t".join(map(str,[acc,taxid,source])))
    seen.add((acc,taxid,source,sourceid))
