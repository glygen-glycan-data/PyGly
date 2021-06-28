#!/bin/env python2

import sys, time, traceback
from collections import defaultdict

from getwiki import GlycanData, Glycan
w = GlycanData()

import findpygly
from pygly.GlycanResource import GlyTouCan
from pygly.GlycanResource import GlyCosmos

def accessions(args):
    if len(args) == 0:
	for it in sys.stdin:
	    yield it.strip()
    else:
	for fn in args:
	    for it in open(fn):
		yield it.strip()

gtc = GlyTouCan(usecache=False)
gco = GlyCosmos(usecache=False)

allgco = set(gco.allaccessions())

# allmotifs = dict()
# for acc,label,redend in gtc.allmotifs():
#     allmotifs[acc] = dict(label=label,redend=redend)


archived = set(map(lambda d: d['accession'],gco.archived()))
print "%d accessions archived."%(len(archived),)


current = set()
for gtcacc in accessions(sys.argv[1:]):
    start = time.time()

    g = w.get(gtcacc)
    newgly = False

    if gtcacc in archived:
	if g:
	    w.delete(gtcacc)
	    print >>sys.stderr, "%s deleted in %.2f sec"%(g.get('accession'),time.time()-start,)
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
    iupacseq = gco.getseq(gtcacc,'iupac_extended')
    if iupacseq != None and "," not in iupacseq:
        g.set_annotation(value=iupacseq,
                         property='IUPAC',
                         source='GlyCosmos',type='Sequence')
    else:
	g.delete_annotations(property='IUPAC',source='GlyCosmos',type='Sequence')

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
    g.set_annotation(value=list(set(gtc.gettaxa(gtcacc))),
                     property='Taxonomy', sourceid=gtcacc, 
                     source='GlyTouCan', type='Taxonomy')

    g.delete_annotations(source='GlyTouCan',type='Publication')
    g.set_annotation(value=list(set(gtc.getrefs(gtcacc))),
                     property='Publication', sourceid=gtcacc, 
                     source='GlyTouCan', type='Publication')

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
	    print >>sys.stderr, "%s created in %.2f sec"%(g.get('accession'),time.time()-start,)
	else:
	    print >>sys.stderr, "%s updated in %.2f sec"%(g.get('accession'),time.time()-start,)
    else:
	print >>sys.stderr, "%s checked in %.2f sec"%(g.get('accession'),time.time()-start,)

    current.add(gtcacc)

