#!/bin/env apython3

import sys, csv, re

# ProteinName   PeptideSequence Glycans Mods    RetentionTime
# NormalizedRetentionTime PrecursorMz     PrecursorCharge ProductMz
# ProductCharge   Annotation      LibraryIntensity        Scans

from collections import defaultdict

scans = defaultdict(lambda: defaultdict(list))
for fn in sys.argv[1:]:
  base = fn.rsplit('.',4)[0]
  for row in csv.DictReader(open(fn),dialect='excel-tab'):
    gpkey = tuple(map(row.get,("PeptideSequence","Glycans","Mods","PrecursorCharge")))
    for sc in row.get('Scans').split(';'):
        scdata = dict(zip(("scan","mode","rt","precscan","precmz","score"),re.split(r'[(),]',sc)))
	pmz = float(row['PrecursorMz'])
	pz = int(row['PrecursorCharge'])
	pm = pmz*pz-1.0078*pz
        scoredata = dict(Spectrum=dict(PrecursorMz="%.3f"%pmz,PrecursorMass="%.3f"%pm))
	scoredata['Glycopeptide'] = dict(Hash=row['GlycopeptideHash'])
        if scdata.get('score').strip():
            for score in scdata['score'].split('/'):
                ss = re.split('[:=]',score)
                scoredata[ss[0]] = dict()
                for i in range(1,len(ss),2):
                    scoredata[ss[0]][ss[i]] = float(ss[i+1])
            scans[(base,scdata['scan'])][gpkey] = scoredata

headers = ["File","Scan","Glycopeptide:Hash","Peptide","Glycan","Mods","Charge","Spectrum:PrecursorMass","GPS:ICScore","Byonic:q-value","MSFragger:E-value","MSFragger:PeptideProphet"]
print("\t".join(headers))
for sc in sorted(scans,key=lambda t: (t[0],int(t[1]))):
    for gpkey in scans[sc]:
        row = dict(File=sc[0],Scan=sc[1])
        row["Peptide"] = gpkey[0]
        row["Glycan"] = gpkey[1]
        row["Mods"] = gpkey[2]
        row["Charge"] = gpkey[3]
        for score in scans[sc][gpkey]:
            for k,v in scans[sc][gpkey][score].items():
                row["%s:%s"%(score,k)] = v
        print("\t".join(map(lambda h: str(row.get(h,"-")),headers)))
