#!/bin/env python2

import sys, time, traceback, hashlib
from collections import defaultdict

import findpygly
from pygly.GlycanResource import GlyTouCanNoCache, GlyTouCan

gtc = GlyTouCanNoCache();

acc2hash = defaultdict(set)

# for acc,inhsh,vwhsh,hgshsh in gtc.query_validacc1():
# if hgshsh != vwhsh:
# 	continue
# acc2hash[acc].add(vwhsh)
# acc2hash[acc].add(inhsh)

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

    if not g.has_annotations(property='GlycoCT',type='Sequence',source='GlyTouCan'):
	if glycan and not g.has_annotations(property='GlycoCT',type='Sequence',source='EdwardsLab'):
	    value = glycan.glycoct()
	    g.set_annotation(property='GlycoCT',type='Sequence',value=value,source='EdwardsLab')
    else:
	g.delete_annotations(source='EdwardsLab',type='Sequence',property='GlycoCT')

    hashes = acc2hash[g.get('accession')]
    for ann in g.annotations(type='Sequence'):
	# if ann.get('property') in ('SequenceHash','IUPAC','GlycoWorkBench','GLYCAM-IUPAC'):
	#     continue
	if ann.get('property') not in ('GlycoCT','WURCS'):
	    continue
	value = ann.get('value')
	if not value:
	    continue
	hashes.add(hashlib.sha256(value).hexdigest().lower())
    if len(hashes) > 0:
	g.set_annotation(property='SequenceHash',type='Sequence',value=sorted(hashes),source='EdwardsLab')

    try:
	if glycan and glycan.fully_determined():
            glycam = glycan.glycam()
	    if glycam and '?' not in glycam:
                g.set_annotation(value=glycam,property="GLYCAM-IUPAC",type='Sequence',source='EdwardsLab')
	    else:
                g.delete_annotations(property="GLYCAM-IUPAC",type='Sequence',source='EdwardsLab')
	else:
            g.delete_annotations(property="GLYCAM-IUPAC",type='Sequence',source='EdwardsLab')
    except (KeyError,AssertionError):
        g.delete_annotations(property="GLYCAM-IUPAC",type='Sequence',source='EdwardsLab')
    except:
        g.delete_annotations(property="GLYCAM-IUPAC",type='Sequence',source='EdwardsLab')
        traceback.print_exc()

    if glycan and not glycan.has_root():
	g.delete_annotations(property="IUPAC",type="Sequence")

    if g.has_annotations(property='IUPAC',type='Sequence',source='GlyCosmos'):
	seq = g.get_annotation_value(property='IUPAC',type='Sequence',source='GlyCosmos')
	if ',' not in seq:
	    # These rules supplied by Sriram Neelamegham <neel@buffalo.edu>
            seq=seq.replace('(','-(')           # sometimes '-' is missing before '('
            seq=seq.replace('--','-')           # sometimes '-' is missing before '('
            seq=seq.replace(')-D',')-?-D')      # sometimes bond type is missing
            seq=seq.replace(')?',')-?')         # sometimes bond type is missing
            seq=seq.replace(')-L',')-?-L')      # sometimes bond type is missing
            seq=seq.replace(']-D',']-?-D')      # sometimes bond type is missing
            seq=seq.replace(']-L',']-?-L')      # sometimes bond type is missing
            seq=seq.replace(']?',']-?')         # sometimes bond type is missing
            seq=seq.replace(')alpha',')-alpha') # sometimes '-' missing before bond
            seq=seq.replace(')beta',')-beta')   # sometimes '-' missing before bond 
	    g.set_annotation(value=seq,property='IUPAC',type='Sequence',source='EdwardsLab')

    g.delete_annotations(source='EdwardsLab',type='MolWt')
    try:
        if glycan:
            mw = glycan.underivitized_molecular_weight()
            g.set_annotation(value=mw,
                             property='UnderivitizedMW',
                             source='EdwardsLab', type='MolWt')
    except (KeyError,ValueError):
	pass
    except:
        traceback.print_exc()

    try:
        if glycan:
            pmw = glycan.permethylated_molecular_weight()
            g.set_annotation(value=pmw,
                             property='PermethylatedMW',
                             source='EdwardsLab', type='MolWt')
    except (KeyError,ValueError):
        pass
    except:
        traceback.print_exc()

    g.delete_annotations(source='EdwardsLab',type='MonosaccharideCount')
    hasmonosaccharides = set()
    try: 
	if glycan:
	    if not glycan.repeated():
                comp = glycan.iupac_composition()
                comp1 = glycan.iupac_composition(floating_substituents=False)
		comp2 = defaultdict(int)
		comp3 = defaultdict(int)
            else:
                comp = glycan.iupac_composition(repeat_times=1)
                comp1 = glycan.iupac_composition(repeat_times=1,floating_substituents=False)
                comp2 = glycan.iupac_composition(repeat_times=2)
                comp3 = glycan.iupac_composition(repeat_times=2,floating_substituents=False)
		for key in comp2:
		    comp2[key] -= comp[key]
		for key in comp3:
		    comp3[key] -= comp1[key]
	else:
	    comp = {}
	    comp1 = {}
	    comp2 = defaultdict(int)
	    comp3 = defaultdict(int)
	for ckey,count in comp.items():
            if count > 0 and not ckey.startswith('_'):
		if ckey=='Count':
		    g.set_annotation(value=(count if comp2[ckey] == 0 else str(count)+"+"),
		         property='MonosaccharideCount',
		         source='EdwardsLab',type='MonosaccharideCount')
		else:
	            g.set_annotation(value=(count if comp2[ckey] == 0 else str(count)+"+"),
		        property=ckey+'Count',
		        source='EdwardsLab',type='MonosaccharideCount')
		    hasmonosaccharides.add(ckey)
	for ckey,count in comp1.items():
	    if count > 0 and not ckey.startswith('_') and ckey != "Count":
		if ckey.endswith('+aldi'):
		    g.set_annotation(value=(count if comp3[ckey] == 0 else str(count)+"+"),
			property=ckey+"Count",
			source='EdwardsLab',type='MonosaccharideCount')
		    hasmonosaccharides.add(ckey)

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
	    g.set_annotation(value=glycan.iupac_redend(),property='ReducingEnd',
                             source='EdwardsLab',type='Structure')
	    g.set_annotation(value=hasmonosaccharides,property='HasMonosaccharide',
                             source='EdwardsLab',type='Structure')
    except:
        traceback.print_exc()

    try:
        glycoctseq = g.get_annotation_value(property='GlycoCT',type='Sequence')
    except LookupError:
        glycoctseq = None
    try:
        wurcsseq = g.get_annotation_value(property='WURCS',type='Sequence')
    except LookupError:
        wurcsseq = None
    repstr = False
    if (wurcsseq and '~' in wurcsseq) or (glycoctseq and 'REP' in glycoctseq):
        repstr = True
    g.set_annotation(value=str(repstr).lower(),property='HasRepeat',source='EdwardsLab',type='Structure')

    if w.put(g):
        print >>sys.stderr, "%s updated in %.2f sec"%(g.get('accession'),time.time()-start,)
    else:
	print >>sys.stderr, "%s checked in %.2f sec"%(g.get('accession'),time.time()-start,)

