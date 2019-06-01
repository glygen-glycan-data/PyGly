#!/bin/env python27

import sys, time, traceback
from collections import defaultdict

from getwiki import GlycanDataWiki, Glycan
w = GlycanDataWiki()

import findpygly
from pygly.GlyTouCan import GlyTouCan

remove_extras = False
if len(sys.argv) >= 2 and sys.argv[1] == "--clean":
    remove_extras = True

def accessions(args):
    if len(args) == 0:
	for it in sys.stdin:
	    yield it.strip()
    else:
	for fn in args:
	    for it in open(fn):
		yield it.strip()

gtc = GlyTouCan(usecache=True)

current = set()
for gtcacc in accessions(sys.argv[1:]):
    start = time.time()

    g = w.get(gtcacc)
    newgly = False

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
    g.set_annotation(value=gtc.getseq(gtcacc,'iupac_extended'),
                     property='IUPAC',
                     source='GlyTouCan',type='Sequence')

    g.delete_annotations(source='GlyTouCan',type='MolWt')
    g.set_annotation(value=gtc.getmass(gtcacc),
                     property='UnderivitizedMW',
                     source='GlyTouCan',type='MolWt')

    g.delete_annotations(source='GlyTouCan',type='MonosaccharideCount')
    g.set_annotation(value=gtc.getmonocount(gtcacc),
		             property='MonosaccharideCount',
		             source='GlyTouCan',type='MonosaccharideCount')

    g.delete_annotations(source='GlyTouCan',type='CrossReference')
    dic = defaultdict(list)
    xref_dic = {'glycosciences_de':'GLYCOSCIENCES.de',
                'pubchem':'PubChem',
                'kegg':'KEGG',
                'unicarbkb':'UniCarbKB',
                'glyconnect':'GlyConnect',
                'glycome-db':'GlycomeDB',
                'cfg':'CFG',
                'pdb':'PDB',
		'bcsdb':'BCSDB',
                'carbbank':'Carbbank(CCSB)'}    
    for xref in gtc.getcrossrefs(gtcacc):
        ref, c = xref.split(":")
	dic[ref].append(c)
    for key in dic:       
        g.set_annotation(value=dic[key],
                         property=xref_dic.get(key,key),
                         source='GlyTouCan',type='CrossReference')
    g.set_annotation(value=gtcacc,property='GlyTouCan',
                     source='GlyTouCan',type='CrossReference')

    g.delete_annotations(source='GlyTouCan',type='Motif')
    g.set_annotation(value=gtc.getmotif(gtcacc),
                     property='Motif',
                     source='GlyTouCan', type='Motif')

    g.delete_annotations(source='GlyTouCan',type='Taxonomy')
    g.set_annotation(value=list(set(gtc.gettaxa(gtcacc))),
                     property='Taxonomy',
                     source='GlyTouCan', type='Taxonomy')

    g.delete_annotations(source='GlyTouCan',type='Publication')
    g.set_annotation(value=list(set(gtc.getrefs(gtcacc))),
                     property='Publication',
                     source='GlyTouCan', type='Publication')

    g.delete_annotations(source='GlyTouCan',type='Subsumption')
    topo = gtc.gettopo(gtcacc)
    if topo:
        g.set_annotation(value=topo,
                         property='Topology',
                         source='GlyTouCan', type='Subsumption')
    comp = gtc.getcomp(gtcacc)
    if comp:
        g.set_annotation(value=comp,
                         property='Composition',
                         source='GlyTouCan', type='Subsumption')
    basecomp = gtc.getbasecomp(gtcacc)
    if basecomp:
        g.set_annotation(value=basecomp,
                         property='BaseComposition',
                         source='GlyTouCan', type='Subsumption')

    if gtcacc == topo:
        g.set_annotation(value='Topology',
                         property='SubsumptionLevel',
                         source='GlyTouCan', type='Subsumption')
    elif gtcacc == comp:
        g.set_annotation(value='Composition',
                         property='SubsumptionLevel',
                         source='GlyTouCan', type='Subsumption')
    elif gtcacc == basecomp:
        g.set_annotation(value='BaseComposition',
                         property='SubsumptionLevel',
                         source='GlyTouCan', type='Subsumption')
    else:
        g.set_annotation(value='Saccharide',
                         property='SubsumptionLevel',
                         source='GlyTouCan', type='Subsumption')
    
    if w.put(g):
	if newgly:
	    print >>sys.stderr, "%s created in %.2f sec"%(g.get('accession'),time.time()-start,)
	else:
	    print >>sys.stderr, "%s updated in %.2f sec"%(g.get('accession'),time.time()-start,)
    else:
	print >>sys.stderr, "%s checked in %.2f sec"%(g.get('accession'),time.time()-start,)

    current.add(gtcacc)

if remove_extras:
    for g in w.iterglycan():
        if g.get('accession') not in current:
            print >>sys.stderr, "Deleting:",g.get('pagename')
            g.delete(w)
