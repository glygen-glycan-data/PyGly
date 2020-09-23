#!/bin/env apython

import sys, os, os.path, json, csv, array
from collections import defaultdict
import getwiki
from analysis.fdr import CombinedAnalysisFDR
from analysis.regression import SimpleLinearRegression as SLR

from optparse import OptionParser

parser = OptionParser()
parser.add_option("--transitions",type='string',dest='transitions',default=None,
                  help="Glycopeptide transitions analyzed. Required.")
parser.add_option("--chromatograms",type='string',dest='chromatograms',default=None,
                  help="Glycopeptide chromatograms output by OpenSWATH. Required.")
parser.add_option("--results",type='string',dest='results',default=None,
                  help="Glycopeptide results table output by OpenSWATH. Required.")
parser.add_option("--ndecoys",type="int",dest='ndecoys',default=1,
                  help="Number of decoy transitions per target transition. Default: 1.")
parser.add_option("--thresh",type="float",dest='thresh',default=1,
                  help="FDR treshold (in %) for results filtering. Default: 1%.")
parser.add_option("--outdir",type='string',dest='outdir',default=None,
                  help="Output directory for JSON chromatograms. Required.")

opts,args = parser.parse_args()                                                                                              

if not opts.transitions:
    parser.error("Glycopeptide transitions file missing.")
if not opts.chromatograms:
    parser.error("Glycopeptide chromatograms from OpenSWATH missing.")
if not opts.results:
    parser.error("Glycopeptide results table from OpenSWATH missing.")
if not opts.outdir:
    parser.error("Output directory for JSON chromatograms missing.")

tgkey = "transition_group_id"
scorekey = "main_var_xx_swath_prelim_score"

fdr = CombinedAnalysisFDR(rankingkey=scorekey,ndecoys=opts.ndecoys)
fdr.add_ids(csv.DictReader(open(opts.results),dialect="excel-tab"))
# for sc,qv,ta,de in fdr.allqvalues():
#     print sc,qv,ta,de
scorethreshold = fdr.score(opts.thresh/100.0)
targets = fdr.targets(scorethreshold)
print "Using",targets,"results for score threshold:",scorethreshold,"for",opts.thresh,"% FDR"

alltr = dict()
pepgrp = defaultdict(list)

last_transition_group_id = None
for row in sorted(csv.DictReader(open(opts.results),dialect="excel-tab"),
                  key=lambda r: (r[tgkey],-float(r[scorekey]))):
    transition_group_id = row[tgkey]
    if transition_group_id == last_transition_group_id:
        continue
    last_transition_group_id = transition_group_id
    decoy = int(row['decoy'])
    if decoy == 1:
	continue
    score = float(row[scorekey])
    if score < scorethreshold:
        continue
    apex = map(float,row['aggr_Peak_Apex'].split(';'))
    area = map(float,row['aggr_Peak_Area'].split(';'))
    trs = row['aggr_Fragment_Annotation'].split(';')
    rt = float(row['RT'])/60.0
    assay_rt = float(row['assay_rt'])/60.0
    for tr,ap,ar in zip(trs,apex,area):
        alltr[tr] = dict(area=ar,apex=ap,tg=transition_group_id,exprt=rt,transnrt=assay_rt)

for row in csv.DictReader(open(opts.transitions),dialect="excel-tab"):
    tr = row["transition_name"]
    pepid = row["peptide_id"]
    pmz = float(row["PrecursorMz"])
    pz = int(row["PrecursorCharge"])
    prmz = float(row["ProductMz"])
    prz = int(row["ProductCharge"])
    nrt = float(row["Tr_recalibrated"])
    libint = int(row["LibraryIntensity"])
    if tr in alltr:
        alltr[tr].update(dict(pepid=pepid,pmz=pmz,pz=pz,prmz=prmz,prz=prz,nrt=nrt,libint=libint))
        pepgrp[(pepid,pz)].append(tr)

nrtpairs = []
exprt=[]
for tr in alltr:
    x = alltr[tr]['transnrt']
    y = alltr[tr]['nrt']
    x1 = alltr[tr]['exprt']
    nrtpairs.append((x,y))
    exprt.append(x1)
regress = SLR()
params = regress.fit(nrtpairs)
print "Exp. RT = %f * NRT %+f"%(1/params['slope'],-params['intercept']/params['slope'])
for tr in alltr:
    alltr[tr]['expnrt'] = regress.y(params,alltr[tr]['exprt'])

try:
    import xml.etree.cElementTree as ET
except ImportError:
    try:
        import xml.etree.ElementTree as ET
    except:
        import elementtree.ElementTree as ET

import PyMSNumpress, zlib, base64

def maketag(*tags):
    return "/".join(['{http://psi.hupo.org/ms/mzml}'+t for t in tags])

def msnumpress_decode(fn,s):
    result = []
    fn(s,result)
    return result

context = []
for event,ele in ET.iterparse(open(opts.chromatograms),('start','end')):
    if event == "start":
        context.append(ele.tag)
    elif event == "end":
        del context[-1]

    if event == "end" and ele.tag == maketag("chromatogram"):
        if ele.attrib['id'] not in alltr:
            continue
        rt = None; intensity = None
        bdal = ele.find(maketag('binaryDataArrayList'))
        for bda in bdal.findall(maketag('binaryDataArray')):
            convert = None
            axis = None
            for cvp in bda.findall(maketag('cvParam')):
                if cvp.attrib['accession'] == "MS:1000595":
                    axis = 'rt'
                if cvp.attrib['accession'] == "MS:1000515":
                    axis = 'int'
                if cvp.attrib['accession'] == "MS:1002746":
                    # "MS-Numpress linear prediction compression followed by zlib compression"
                    convert = (lambda s: msnumpress_decode(PyMSNumpress.decodeLinear,
                                                           array.array('B',
                                                                       zlib.decompress(base64.b64decode(s)))))
                if cvp.attrib['accession'] == "MS:1002748":
                    # "MS-Numpress short logged float compression followed by zlib compression"
                    convert = (lambda s: msnumpress_decode(PyMSNumpress.decodeSlof,
                                                           array.array('B',
                                                                       zlib.decompress(base64.b64decode(s)))))

            assert convert != None and axis in ('rt','int')
            values = convert(bda.find(maketag('binary')).text)
            if axis == "rt":
                values = map(lambda x: x/60.0, values)
            alltr[ele.attrib['id']][axis] = values

if not os.path.isdir(opts.outdir):
    os.makedirs(opts.outdir)

for i,(pepid,charge) in enumerate(pepgrp):
    fn = '%s.%s.50.json'%(pepid,charge)
    data = {}
    data["trlib_z1"] = charge
    data["trlib_pepid"] = pepid
    # data["trlib_mz1"] = None
    # data["rt"] = None
    # data["nrt"] = None
    data["xaxis"]='rt'
    data["series"] = []
    maxarea = max([alltr[trid]['area'] for trid in pepgrp[(pepid,charge)]])
    for trid in pepgrp[(pepid,charge)]:
        tr = alltr[trid]
        trdata = dict()
        trdata["trlib_relint"] = tr['libint']
        trdata["obs_relint"] = int(round(1000.0*tr['area']/maxarea,0))
        trdata["trlib_id"] = ".".join(map(str,[pepid,charge,trid.split("_")[-1]]))
        # trdata["name"] = None
        trdata["yaxis"] = {"label": "int"}
        trdata["pairs"] = zip(tr["rt"],tr["int"])
        if "nrt" not in data:
            data["nrt"] = tr["expnrt"]
        if "rt" not in data:
            data["rt"] = tr["exprt"]
        if "trlib_mz1" not in data:
            data["trlib_mz1"] = tr["pmz"]
        trdata["name"] = trid
        data["series"].append(trdata)
    # print i,pepid,charge
    wh = open(os.path.join(opts.outdir,fn),'w')
    wh.write(json.dumps(data))
    wh.close()

print "%d peptide groups extracted"%(len(pepgrp),)
