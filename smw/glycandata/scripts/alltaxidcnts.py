#!/bin/env python3.12

import sys, time, traceback, csv
from collections import defaultdict

import findpygly
from pygly.GlycanResource import GlyTouCan
from pygly.GlycanResource import GlyCosmos
from pygly.GlycanResource import GlyConnect
from pygly.Taxonomy import NCBITaxonomy

gtc = GlyTouCan(verbose=False,usecache=False)
gco = GlyCosmos(verbose=False,usecache=False)
gcn = GlyConnect(verbose=False)
tax = NCBITaxonomy()

archived = set(map(lambda d: d['accession'],gco.archived()))
# print("%d accessions archived."%(len(archived),))


def gentaxa():
    for r in csv.DictReader(open(sys.argv[1]),dialect='excel-tab'):
        yield r['GlyTouCanAccession'],r['TaxID'],'GlyGen'
    for gtcacc,taxid in gtc.alltaxa():
        yield gtcacc,taxid,'GlyTouCan'
    for gtcacc,taxid in gco.alltaxa():
        yield gtcacc,taxid,'GlyCosmos'
    for gcnid,gtcacc,taxid in gcn.alltaxa():
        yield gtcacc,taxid,'GlyConnect'

freq = defaultdict(set)
alltaxid = set()

for gtcacc,taxid,src in gentaxa():
    if gtcacc in archived:
        continue
    try:
        alltaxid.add(int(taxid))
    except ValueError:
        continue
    freq[(src,'as annotated',int(taxid))].add(gtcacc)
    ranks = tax.get_ranked_ancestors(taxid)
    for r in ('species','genus','family'):
        if r in ranks:
            alltaxid.add(ranks[r])
            freq[(src,r,ranks[r])].add(gtcacc)

print("\t".join(["NCBI TaxID","Scientific Name","Common Name","Rank","Aggr","GlyTouCan","GlyCosmos","GlyConnect","GlyGen","Total"]));
for aggr in ("as annotated","species","genus"):
    for taxid in sorted(alltaxid):
        allgtc = set()
        for src in ("GlyTouCan","GlyCosmos","GlyConnect","GlyGen"):
            allgtc.update(freq[(src,aggr,taxid)])
        if len(allgtc) == 0:
            continue
        print("\t".join(map(str,[taxid,
                                 tax.get_scientific_name(taxid),
                                 tax.get_common_name(taxid,""),
                                 tax.get_rank(taxid,""),
                                 aggr,
                                 len(freq[("GlyTouCan",aggr,taxid)]),
                                 len(freq[("GlyCosmos",aggr,taxid)]),
                                 len(freq[("GlyConnect",aggr,taxid)]),
                                 len(freq[("GlyGen",aggr,taxid)]),
                                 len(allgtc)])))
