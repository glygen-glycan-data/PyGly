#!/bin/env python3.12

import sys, time, traceback
from collections import defaultdict

from getwiki import GlycanData, Glycan
w = GlycanData()

import findpygly
from pygly.GlycanResource import GlyTouCan, GlyCosmos
from pygly.GlycanResource import GlyProtDBSourceFile, GlycomeAtlasSourceFile

def accessions(args):
    if len(args) == 0:
        for it in sys.stdin:
            yield it.strip()
    else:
        for fn in args:
            for it in open(fn):
                yield it.strip()

gtc = GlyTouCan(verbose=False,prefetch=True,usecache=False)
gco = GlyCosmos(verbose=False,prefetch=True,usecache=False)

allgco = set(gco.allaccessions())

# allmotifs = dict()
# for acc,label,redend in gtc.allmotifs():
#     allmotifs[acc] = dict(label=label,redend=redend)

archived = set(map(lambda d: d['accession'],gco.archived()))
print("%d accessions archived."%(len(archived),),file=sys.stderr)

# get GlyCosmos taxa entries from GlyGen Source Files

gpdb = GlyProtDBSourceFile()
glyat = GlycomeAtlasSourceFile()
glygentaxa = defaultdict(set)
for clsinst in (gpdb,glyat):
  for sourceacc,gtcacc,taxid,source,sourceid in clsinst.alltaxa():
    if source != "GlyCosmos":
        continue
    glygentaxa[gtcacc].add(taxid)

current = set()
for gtcacc in accessions(sys.argv[1:]):
    start = time.time()

    g = w.get(gtcacc)
    newgly = False

    if gtcacc in archived:
        if g:
            w.delete(gtcacc)
            print("%s deleted in %.2f sec"%(g.get('accession'),time.time()-start,),file=sys.stderr)
        continue

    if not g:
        newgly = True
        g = Glycan(accession=gtcacc)

    g.delete('wurcs')
    g.delete('glycoct')
    g.delete('iupac')

    g.delete_annotations(source='GlyTouCan',type='Sequence')
    g.set_annotation(value=gtc.getseq(gtcacc,'wurcs'),
                     property='WURCS',
                     source='GlyTouCan',type='Sequence')
    g.set_annotation(value=gtc.getseq(gtcacc,'glycoct'),
                     property='GlycoCT',
                     source='GlyTouCan',type='Sequence')

    g.delete_annotations(property='IUPAC',source='GlyTouCan',type='Sequence')
    g.delete_annotations(property='IUPAC',source='GlyCosmos',type='Sequence')

    iupacseq = gco.getseq(gtcacc,'iupac_extended')
    if iupacseq != None and "," not in iupacseq and "'" not in iupacseq:
        try:
            # iupacseq = iupacseq.toPython()
            # iupacseq1 = iupacseq.replace(u'\u2192','->').replace(u'\u2194','<->').replace(u'\u03b1','alpha').replace(u'\u03b2','beta')
            # iupacseq2 = iupacseq.replace(u'\u2192','->').replace(u'\u2194','<->').replace(u'\u03b1','a').replace(u'\u03b2','b')
            pass
        except AttributeError:
            # print(gtcacc,repr(iupacseq))
            # raise 
            pass
        g.set_annotation(value=iupacseq, property='IUPAC-Extended', source='GlyCosmos',type='Sequence')
    else:
        g.delete_annotations(property='IUPAC-Extended',source='GlyCosmos',type='Sequence')

    iupacseq = gco.getseq(gtcacc,'iupac_condensed')
    if iupacseq != None and "," not in iupacseq and "'" not in iupacseq:
        g.set_annotation(value=iupacseq, property='IUPAC-Condensed', source='GlyCosmos',type='Sequence')
    else:
        g.delete_annotations(property='IUPAC-Condensed',source='GlyCosmos',type='Sequence')

    g.delete_annotations(source='GlyTouCan',type='MolWt')
    g.set_annotation(value=gtc.getmass(gtcacc),
                     property='UnderivitizedMW',
                     source='GlyTouCan',type='MolWt')

    g.delete_annotations(source='GlyTouCan',type='MonosaccharideCount')
    g.set_annotation(value=gtc.getmonocount(gtcacc),
                             property='MonosaccharideCount',
                             source='GlyTouCan',type='MonosaccharideCount')

    xref_dic = {# 'glycosciences_de':'GLYCOSCIENCES.de',
                # 'pubchem':'PubChem',
                # 'kegg':'KEGG',
                'kegg_glycan':'KEGG',
                # 'unicarbkb':'UniCarbKB',
                'unicarb-db':'UniCarb-DB',
                'glyconnect':'GlyConnectStructure',
                'glyconnect-comp':'GlyConnectComposition',
                # 'glycome-db':'GlycomeDB',
                # 'cfg':'CFG',
                # 'pdb':'PDB',
                'bcsdb':'BCSDB',
                'matrixdb':'MatrixDB',
                'glycoepitope':'GlycoEpitope',
                # 'carbbank':'Carbbank(CCSB)',
               }    
    for prop in xref_dic.values():
        g.delete_annotations(source='GlyTouCan',property=prop,type='CrossReference')
    dic = defaultdict(list)
    for ref, c in gtc.getcrossrefs(gtcacc):
        # ref, c = xref.split(":")
        dic[ref].append(c)
    for key in dic:       
        if key not in xref_dic:
            continue
        g.set_annotation(value=dic[key],
                         property=xref_dic[key],
                         source='GlyTouCan',type='CrossReference')
    g.set_annotation(value=gtcacc,property='GlyTouCan',
                     source='GlyTouCan',type='CrossReference')
    if gtcacc in allgco:
        g.set_annotation(value=gtcacc,property='GlyCosmos',
                         source='GlyCosmos',type='CrossReference')

    # g.delete_annotations(source='GlyTouCan',type='Motif')
    # g.set_annotation(value=gtc.getmotif(gtcacc),
    #                  property='Motif',
    #                  source='GlyTouCan', type='Motif')
    # value = map(lambda acc: ":".join([acc,allmotifs[acc]['label'],allmotifs[acc]['redend']]),gtc.getmotif(gtcacc))
    # g.set_annotation(value=value,
    #                  property='NamedMotif',
    #                  source='GlyTouCan', type='Motif')

    g.delete_annotations(source='GlyTouCan',type='Taxonomy')
    # g.set_annotation(value=list(set(gtc.gettaxa(gtcacc))),
    #                  property='Taxonomy', sourceid=gtcacc, 
    #                  source='GlyTouCan', type='Taxonomy')

    g.delete_annotations(source='GlyCosmos',type='Taxonomy')
    taxmap = {'11103': '3052230'}
    taxids = set(gco.gettaxa(gtcacc))
    taxids.update(glygentaxa[gtcacc])
    taxids = set(map(lambda t: taxmap.get(t,t),taxids))
    g.set_annotation(value=list(sorted(taxids)), property='Taxonomy', sourceid=gtcacc, 
                     source='GlyCosmos', type='Taxonomy')

    # g.delete_annotations(source='GlyTouCan',type='Publication')
    # g.set_annotation(value=list(set(gtc.getrefs(gtcacc))),
    #                  property='Publication', sourceid=gtcacc, 
    #                  source='GlyTouCan', type='Publication')

    # get this stuff from GNOme, now...
    g.delete_annotations(source='GlyTouCan',type='Subsumption')

    # topo = gtc.gettopo(gtcacc)
    # if topo:
    #    g.set_annotation(value=topo,
    #                     property='Topology',
    #                     source='GlyTouCan', type='Subsumption')
    # comp = gtc.getcomp(gtcacc)
    # if comp:
    #     g.set_annotation(value=comp,
    #                      property='Composition',
    #                      source='GlyTouCan', type='Subsumption')
    # basecomp = gtc.getbasecomp(gtcacc)
    # if basecomp:
    #     g.set_annotation(value=basecomp,
    #                      property='BaseComposition',
    #                      source='GlyTouCan', type='Subsumption')

    # if gtcacc == topo:
    #     g.set_annotation(value='Topology',
    #                      property='SubsumptionLevel',
    #                      source='GlyTouCan', type='Subsumption')
    # elif gtcacc == comp:
    #     g.set_annotation(value='Composition',
    #                      property='SubsumptionLevel',
    #                      source='GlyTouCan', type='Subsumption')
    # elif gtcacc == basecomp:
    #     g.set_annotation(value='BaseComposition',
    #                      property='SubsumptionLevel',
    #                      source='GlyTouCan', type='Subsumption')
    # else:
    #     g.set_annotation(value='Saccharide',
    #                      property='SubsumptionLevel',
    #                      source='GlyTouCan', type='Subsumption')
    
    if w.put(g):
        if newgly:
            print("%s created in %.2f sec"%(g.get('accession'),time.time()-start,),file=sys.stderr)
        else:
            print("%s updated in %.2f sec"%(g.get('accession'),time.time()-start,),file=sys.stderr)
    else:
        print("%s checked in %.2f sec"%(g.get('accession'),time.time()-start,),file=sys.stderr)

    current.add(gtcacc)

