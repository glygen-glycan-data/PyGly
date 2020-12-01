#!/bin/env python27

import re, sys

import findpygly
from pygly.GlycanResource import GlyGen, GlyCosmos

glygen = GlyGen(usecache=False)
glycosmos = GlyCosmos(usecache=False)

allglygen = set()
for acc in glygen.allglycans():
    allglygen.add(acc)

archived = set()
for item in glycosmos.archived():
    item['accession'] = str(item['accession'])
    if item['accession'] not in allglygen:
	continue
    archived.add(item['accession'])

replace = dict()
for item in glycosmos.replaced():
    item['replace'] = str(item['replace'])
    item['accession'] = str(item['accession'])
    if item['replace'] not in allglygen:
	continue
    if item['accession'] not in allglygen:
	continue
    assert item['replace'] not in replace
    assert item['replace'] in archived
    replace[item['replace']] = item['accession']

print "\t".join(["accession","replacewith"])
for acc in sorted(archived):
    print "\t".join(filter(None,[acc,replace.get(acc,"-")]))
