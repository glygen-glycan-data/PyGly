#!/bin/env python27

import sys, time, traceback

from getwiki import GlycanData
w = GlycanData()

def accessions():
    if len(sys.argv) > 1:
	for arg in sys.argv[1:]:
	    g = w.get(arg.strip())
	    if g:
		yield g
    else:
	for g in w.iterglycan():
	    yield g

for g in accessions():
    start = time.time()

    glycan = g.getGlycan()

    if not glycan:
	continue

    if not g.has_annotations(property='GlycoCT',type='Sequence',source='GlyTouCan'):
	if not g.has_annotations(property='GlycoCT',type='Sequence',source='EdwardsLab'):
	    value = glycan.glycoct()
	    g.set_annotation(property='GlycoCT',type='Sequence',value=value,source='EdwardsLab')
    else:
	g.delete_annotations(source='EdwardsLab',type='Sequence',property='GlycoCT')

    try:
	if glycan.fully_determined():
            glycam = glycan.glycam()
	    if glycam and '?' not in glycam:
                g.set_annotation(value=glycam,property="GLYCAM-IUPAC",type='Sequence',source='EdwardsLab')
	    else:
                g.delete_annotations(property="GLYCAM-IUPAC",type='Sequence',source='EdwardsLab')
	else:
            g.delete_annotations(property="GLYCAM-IUPAC",type='Sequence',source='EdwardsLab')
    except KeyError:
        g.delete_annotations(property="GLYCAM-IUPAC",type='Sequence',source='EdwardsLab')
    except:
        g.delete_annotations(property="GLYCAM-IUPAC",type='Sequence',source='EdwardsLab')
        traceback.print_exc()

    g.delete_annotations(source='EdwardsLab',type='MolWt')
    try:
        if glycan:
            mw = glycan.underivitized_molecular_weight()
            g.set_annotation(value=mw,
                             property='UnderivitizedMW',
                             source='EdwardsLab', type='MolWt')
    except KeyError:
	pass
    except:
        traceback.print_exc()

    try:
        if glycan:
            pmw = glycan.permethylated_molecular_weight()
            g.set_annotation(value=pmw,
                             property='PermethylatedMW',
                             source='EdwardsLab', type='MolWt')
    except KeyError:
        pass
    except:
        traceback.print_exc()

    g.delete_annotations(source='EdwardsLab',type='MonosaccharideCount')
    try: 
	if glycan:
            comp = glycan.iupac_composition()
	else:
	    comp = {}
	for ckey,count in comp.items():
            if count > 0 and not ckey.startswith('_'):
		if ckey=='Count':
		    g.set_annotation(value=count,
		         property='MonosaccharideCount',
		         source='EdwardsLab',type='MonosaccharideCount')
		else:
	            g.set_annotation(value=count,
		        property=ckey+'Count',
		        source='EdwardsLab',type='MonosaccharideCount')
    except KeyError:
        pass
    except:
        traceback.print_exc()

    g.delete_annotations(source='EdwardsLab',type='Structure')
    try:
	if glycan:
	    g.set_annotation(value='true' if glycan.fully_determined() else 'false',
                             property='FullyDetermined',
                             source='EdwardsLab',type='Structure')
	    g.set_annotation(value='true' if (glycan.undetermined() and glycan.has_root()) else 'false',
                             property='UndeterminedTopology',
                             source='EdwardsLab',type='Structure')
	    g.set_annotation(value='true' if (not glycan.has_root()) else 'false',
                             property='Composition',
                             source='EdwardsLab',type='Structure')
    except:
        traceback.print_exc()

    if w.put(g):
        print >>sys.stderr, "%s updated in %.2f sec"%(g.get('accession'),time.time()-start,)
    else:
	print >>sys.stderr, "%s checked in %.2f sec"%(g.get('accession'),time.time()-start,)

