#!/bin/env python27

import sys,traceback,re
from collections import defaultdict

from getwiki import GlycoMotifWiki, GlyGenMotif
w = GlycoMotifWiki()

motif2gd = defaultdict(set)

from dataset import CSVFileTable
for r in CSVFileTable(sys.argv[1]):
    entry = r['term (main_entry)'].strip()
    xrefs = r['term_xref']
    allmid = set()
    for xr in filter(None,xrefs.split('|')):
	try:
            src,mid = xr.split(':',1)
        except ValueError:
            continue
	if src.lower() == "glycomotif":
            if not re.search(r'^GGM\.\d{6}$',mid):
		print "Bad motif id: %s (%s)"%(mid,entry)
		continue
            motif2gd[mid].add(entry)

for mid in w.site.allpages(prefix='GGM.',generator=False):
    # print mid
    m = w.get(mid)
    entries = motif2gd[mid]
    # print entries
    dbxrefs = m.get('dbxref',set())
    # print dbxrefs
    newdbxrefs = set()
    for acc in entries:
        newdbxrefs.add(('GlycanDictionary',acc))
    for db,acc in dbxrefs:
        if db != 'GlycanDictionary':
            newdbxrefs.add((db,acc))
    # print newdbxrefs
    if len(newdbxrefs) == 0:
        m.delete('dbxref') 
    else:
        m.set('dbxref',newdbxrefs)
    if w.put(m):
        print m.get('id')
