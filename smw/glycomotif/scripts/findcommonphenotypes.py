#!/bin/env python3

import sys
import csv
import gzip
import math
from poibin import PoiBin

from collections import defaultdict

from getwiki import GlycoMotifWiki

motif2ggm = defaultdict()

mappedalign = {
  'Substructure': 'substr',
  'Whole-Glycan': 'whole',
  'Core': 'core',
  'Nonreducing-End': 'nonred'
}

w = GlycoMotifWiki()
for m in w.itermotif(collection="GGM"):
    mgtc = m.get('glytoucan')
    align = m.get('alignment')[0]
    malign = mappedalign[align]
    motif2ggm[(mgtc,malign)] = m.get("id")

motif2enz = defaultdict(lambda: defaultdict(set))
motif2res = defaultdict(set)
motifs = set()
allenz = set()
for row in csv.DictReader(gzip.open(sys.argv[1],'rt'),dialect='excel-tab'):
    m = row['Motif']
    at = row['AlignmentType']
    if (m,at) not in motif2ggm:
        continue
    m = motif2ggm[(m,at)]
    res = row['MotifResidue']
    if res == "-":
        continue
    motif2res[m].add(res)
    enz = row['HumanEnzyme']
    if enz == "-":
        continue
    enz = enz.split(',')
    motif2enz[m][res].update(enz)
    allenz.update(enz)

enz2pheno = defaultdict(set)
pheno2enz = defaultdict(set)
for row in csv.DictReader(open(sys.argv[2],'rt'),dialect='excel-tab'):
    gene = row['gene_symbol']
    hpo_name = row['hpo_name']
    if gene not in allenz:
        continue
    enz2pheno[gene].add(hpo_name)
    pheno2enz[hpo_name].add(gene)

phenoprob = defaultdict(float)
for pheno in pheno2enz:
    phenoprob[pheno] = len(pheno2enz[pheno])/len(allenz)
    for i in range(1,25):
        if phenoprob[pheno] >= 1.0:
            phenoprob[(pheno,i)] = 1.0
        else:
            phenoprob[(pheno,i)] = 1.0-math.exp(i*math.log(1.0-phenoprob[pheno]))

# for pheno in sorted(pheno2enz):
#     print(pheno,end=" ")
#     for i in range(1,25):
#         print("%.3f"%(phenoprob[(pheno,i)],),end=" ")
#     print()

motif2pheno = defaultdict(lambda: defaultdict(set))
for m in motif2enz:
  for res in motif2enz[m]:
    for enz in motif2enz[m][res]:
      motif2pheno[m][res].update(enz2pheno[enz])

motif2phenocnt = defaultdict(lambda: defaultdict(int))
# motif2phenoresset = defaultdict(lambda: defaultdict(lambda: defaultdict(set)))
for m in motif2pheno:
  for res in motif2pheno[m]:
    for pheno in motif2pheno[m][res]:
      motif2phenocnt[m][pheno] += 1

for m in motif2phenocnt:
  for ph in sorted(motif2phenocnt[m],key=lambda ph: -motif2phenocnt[m][ph]):
    cnt = motif2phenocnt[m][ph]
    probs = []
    for res in motif2enz[m]:
      probs.append(phenoprob[(ph,len(motif2enz[m][res]))])
    pb = PoiBin(probs)
    pr = pb.pval(cnt)
    pr = max(pr,1e-100)
    sc = -10*math.log10(pr)
    if sc < 30:
        continue
    print("\t".join(map(str,[m,cnt,len(motif2enz[m]),round(100*cnt/len(motif2enz[m]),3),pr,-10*math.log10(pr),ph])))

sys.exit(0)

for m in motif2pheno:
  for pheno in motif2phenocnt[m]:
    if motif2phenocnt[m][pheno] > 1:
      print(round(100*motif2phenocnt[m][pheno]/len(motif2res[m]),2),motif2phenocnt[m][pheno],len(motif2res[m])," ".join(m),pheno)
      for res in motif2res[m]:
        if len(motif2phenoresset[m][pheno][res]) > 0: 
          print(" ",res,", ".join(sorted(motif2phenoresset[m][pheno][res])))
