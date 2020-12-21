#!/bin/env python27

import sys, re, copy
from getwiki import GPTWiki
import numpy as np
from collections import defaultdict

from optparse import OptionParser

parser = OptionParser()                                                                                   
parser.add_option("--sample",type='string',dest='sample',default=None,help="Sample identifier")                          
parser.add_option("--method",type='string',dest='method',default=None,help="Method identifier")                          
parser.add_option("--acqtype",type='choice',choices=["DDA","DIA"],dest='acqtype',default=None,help="Acquisition type")
parser.add_option("--inst",type='choice',choices=["Orbitrap","TripleTOF"],dest='inst',default=None,help="Instrument")
parser.add_option("--maxtgs",type="int",dest="maxtgs",default=None,help="Max. transition groups")
parser.add_option("--maxtrans",type="int",dest="maxtrans",default=6,help="Max. transitions per transition group")
parser.add_option("--decoys",type="string",dest="decoys",default="37,41,47,53,59",help="Decoy offset(s), comma separated. Default: 37,41,47,53,59.")

opts,args = parser.parse_args()

if not opts.sample and not opts.method and not opts.acqtype and not opts.inst:
    parser.error("Please set sample, method, instrument, and/or acquitision type")

decoys = map(float,opts.decoys.split(','))

headers = """
peptide_id
transition_group_id
transition_id
transition_group_name
transition_name
PrecursorMz
PrecursorCharge
ProductMz
ProductCharge
Tr_recalibrated
LibraryIntensity
PeptideSequence
ProteinName
FullUniModPeptideName
FragmentType
FragmentSeriesNumber
decoy
""".split()

def label2nmono(lab):
    m = re.search(r'\[(.*)\]',lab)
    kvpairs = re.split(r'([A-Z])',m.group(1))
    nmono = 0
    for i in range(1,len(kvpairs),2):
	if kvpairs[i+1] == "":
	    nmono += 1
	else:
	    nmono += int(kvpairs[i+1])
    return nmono

def label2series(lab):
    return lab[0].lower()

w = GPTWiki()

peps = defaultdict(lambda: defaultdict(dict))

for i,tg in enumerate(w.itertgs(sample=opts.sample,acqtype=opts.acqtype,method=opts.method,inst=opts.inst)):
    pep = w.get(tg.get('peptide'))
    pepid = pep.get('id')
    if pep.get('nrt') == None:
	continue
    z1 = tg.get('z1')
    ntrans = len(tg.get('transitions',[]))
    if ntrans == 0:
	continue
    for trid,trint in tg.get('transitions',[]):
	tr = w.get(trid)
	if trid not in peps[(pepid,z1)]:
            peps[(pepid,z1)][trid]['pepname'] = pep.get('name')
            peps[(pepid,z1)][trid]['nrt'] = pep.get('nrt')
            peps[(pepid,z1)][trid]['seq'] = pep.get('sequence')
            peps[(pepid,z1)][trid]['label'] = tr.get('label')
            peps[(pepid,z1)][trid]['nmono'] = label2nmono(tr.get('label'))
            peps[(pepid,z1)][trid]['series'] = label2series(tr.get('label'))
            peps[(pepid,z1)][trid]['mz2'] = tr.get('mz2')
            peps[(pepid,z1)][trid]['z2'] = tr.get('z2')
            peps[(pepid,z1)][trid]['mz1'] = tr.get('mz1')
            peps[(pepid,z1)][trid]['z1'] = tr.get('z1')
            peps[(pepid,z1)][trid]['int'] = []
        peps[(pepid,z1)][trid]['int'].append(trint)
    if opts.maxtgs != None and i >= (opts.maxtgs-1):
	break

for chpep in peps:
    for trid in peps[chpep]:
	peps[chpep][trid]['int'] = np.median(peps[chpep][trid]['int'])

print "\t".join(headers)
for chpep in peps:
    bpi = None
    sint = sorted(map(lambda trid: peps[chpep][trid]['int'],peps[chpep]),reverse=True)
    if (len(sint) < opts.maxtrans) or (sint[opts.maxtrans-1]/sint[0] < 0.1):
	continue
    for i,trid in enumerate(sorted(peps[chpep],key=lambda trid: peps[chpep][trid]['int'], reverse=True)):
	if i >= opts.maxtrans:
	    break
	tr = peps[chpep][trid]
	if bpi == None:
	    bpi = tr['int']
	relint = int(round(1000.0*tr['int']/bpi,0))
	row = {}
	row['peptide_id'] = chpep[0]
	row['transition_group_id'] = "%s/%d"%chpep
	row['transition_group_name'] = "%s/%d"%(tr['pepname'],chpep[1])
	row['transition_id'] = trid
	# row['transition_name'] = "%s/%d"%(tr['label'],tr['z2'])
	row['transition_name'] = "%s/%s/%s"%(row['transition_group_id'],trid,tr['label'])
	row['PrecursorMz'] = tr['mz1']
	row['PrecursorCharge'] = tr['z1']
	row['ProductMz'] = tr['mz2']
	row['ProductCharge'] = tr['z2']
	row['Tr_recalibrated'] = tr['nrt']
	row['LibraryIntensity'] = relint
	row['FullUniModPeptideName'] = tr['seq']
	row['ProteinName'] = chpep[0]
	row['PeptideSequence'] = tr['seq']
	row['FragmentType'] = tr['series']
	row['FragmentSeriesNumber'] = tr['nmono']
        row['decoy'] = 0
	print "\t".join(map(str,map(row.get,headers)))
	for d in decoys:
	    decoy = copy.copy(row)
	    decoy['PrecursorMz'] += d
	    decoy['ProductMz'] += d
	    decoy['transition_name'] = "decoy%d/"%(int(d),)+decoy['transition_name']
	    decoy['transition_group_id'] = "decoy%d/"%(int(d),)+decoy['transition_group_id']
	    decoy['decoy'] = 1
	    print "\t".join(map(str,map(decoy.get,headers)))
