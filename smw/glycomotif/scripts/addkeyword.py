#!/bin/env python3.12

from getwiki import GlycoMotifWiki
import sys, re

w = GlycoMotifWiki()

for m in w.itermotif(collection="GGM"):
    match = False
    if re.search(sys.argv[2],m.get('accession')):
        match = True
    if not match:
        for name in m.get('name',[]):
            if re.search(sys.argv[2],name):
                match = True
                break
    if not match:
        for keyword in m.get('keyword',set()):
            if keyword in sys.argv[2:]:
                match = True
                break
    if match:
        keywords = m.get("keyword",set())
        keywords.add(sys.argv[1])
        m.set("keyword",keywords)
        if w.put(m):
            print(m.get('id'), sys.argv[1], file=sys.stderr)
