#!/bin/env python27

import sys, traceback
import pandas
import os

from getwiki import GlycoMotifWiki, GlydinMotif
w = GlycoMotifWiki()

from pygly.GlyTouCan import GlyTouCan
gtc = GlyTouCan()

fpath = os.path.dirname(os.path.abspath(sys.argv[0])) + "/../data/epitopes.xlsx"
cont_excel = pandas.read_excel(fpath, sheet_name="Sheet1")
cIndex = list(cont_excel.columns)

row = list(cont_excel[cIndex[0]])
glycoct = list(cont_excel[cIndex[1]])

current = set()
for i in range(len(row)):
    rowNum = row[i]
    seq = glycoct[i].strip().replace("\n\n","\n")
    
    try:
        glytoucan, isnew = gtc.register(seq)
    except:
        # traceback.print_exc()
        continue
        
    accession = "r%06d" % int(rowNum)
    motif = GlydinMotif(accession=accession, name=None, glytoucan=glytoucan, redend=None, aglycon=None)
    if w.update(motif):
        print accession
    current.add(accession)

for m in w.itermotif(collection=GlydinMotif):
    if m.get('accession') not in current:
        print "Deleting:", m.get('pagename')
        w.delete(m.get('pagename'))