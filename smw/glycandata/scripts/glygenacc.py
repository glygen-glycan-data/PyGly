#!/bin/env python27

import sys, os, urllib, json, csv
from collections import defaultdict

def glygen_glycan(target,query,api="api"):
    url = 'https://' + api + '.glygen.org/glycan/' + target + '?' + urllib.urlencode(dict(query=json.dumps(query)))
    return json.loads(urllib.urlopen(url).read())

def method1(beta=False):
    query = dict(mass_type="native")
    result = glygen_glycan('search',query,api=("beta-api" if beta else "api"))
    lid = result['list_id']

    query = dict(id=lid,offset=1,limit=1000000,sort="glytoucan_ac",order="asc")
    result = glygen_glycan('list',query,api=("beta-api" if beta else "api"))

    accs1 = set()
    for i,r in enumerate(result['results']):
        accs1.add(r['glytoucan_ac'])

    return accs1

def method2():
    accs2 = set()
    for org in ('human','mouse','rat'):
        for r in csv.DictReader(urllib.urlopen('https://data.glygen.org/ln2wwwdata/reviewed/%s_glycan_accession_sources.csv'%(org,))):
            accs2.add(r['glytoucan_ac'])
    return accs2

def method3():
    accs3 = set()
    for org in ('human','mouse','rat'):
        for r in csv.DictReader(urllib.urlopen('https://data.glygen.org/ln2wwwdata/reviewed/%s_glycan_properties.csv'%(org,))):
            accs3.add(r['glytoucan_ac'])
    return accs3

def method4():
    accs3 = set()
    for org in ('human','mouse','rat'):
        for r in csv.DictReader(urllib.urlopen('https://data.glygen.org/ln2wwwdata/reviewed/%s_glycan_properties.csv'%(org,))):
	    if r['glycan_mass'] and r['glycan_permass']:
                accs3.add(r['glytoucan_ac'])
    return accs3

def method5():
    return open('../data/glygen_accessions.txt').read().split()

freq = defaultdict(int)
sets = []
sets.append(method1())
sets.append(method1(beta=True))
sets.append(method2())
sets.append(method3())
sets.append(method4())
sets.append(method5())

allaccs = set()
for s in sets:
    allaccs.update(s)

for acc in allaccs:
    bm = 0
    for i,s in enumerate(reversed(sets)):
	if acc in s:
	    bm += 10**i
    freq[bm] += 1

# for bm,f in sorted(freq.items()):
#     print "%0*d\t%d"%(len(sets),bm,freq[bm])

for acc in method1():
    print acc
