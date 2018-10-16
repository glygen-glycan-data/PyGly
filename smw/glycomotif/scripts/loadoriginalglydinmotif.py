#!/bin/env python27

import sys, traceback, re
from dataset import XLSXFileTable
import os

from getwiki import GlycoMotifWiki, GlydinMotif
w = GlycoMotifWiki()

from pygly.GlyTouCan import GlyTouCan
gtc = GlyTouCan()

fpath = sys.argv[1]
cont_excel = XLSXFileTable(fpath)

current = set()
for row in cont_excel:
    rowNum = row['ID']
    seq = row['GlycoCT'].strip()
    seq = re.sub('\n\n+','\n',seq)
    
    try:
        glytoucan, isnew = gtc.register(seq)
    except:
        # traceback.print_exc()
        continue
        
    accession = "R%06d" % int(rowNum)
    motif = GlydinMotif(accession=accession, name=None, glytoucan=glytoucan, redend=None, aglycon=None)
    if w.update(motif):
        print accession
    current.add(accession)

for m in w.itermotif(collection=GlydinMotif):
    if m.get('accession') not in current:
        print "Deleting:", m.get('pagename')
        w.delete(m.get('pagename'))
