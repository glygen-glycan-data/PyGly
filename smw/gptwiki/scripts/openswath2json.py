#!/bin/env apython

import sys, os, os.path, json, csv, array, glob, re, shutil
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
parser.add_option("--fdr",type="float",dest='fdr',default=1,
                  help="FDR treshold (in %) for results filtering. Default: 1%.")
parser.add_option("--score",type="float",dest='score',default=1.5,
                  help="Score treshold for results filtering. Default: 1.5.")
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

scorekey = "main_var_xx_swath_prelim_score"

if  os.path.isdir(opts.outdir):
    shutil.rmtree(opts.outdir)

try:
    fdr = CombinedAnalysisFDR(rankingkey=scorekey,ndecoys=opts.ndecoys)
    fdr.add_ids(opts.results)
    scorethreshold = max(opts.score,fdr.score(opts.fdr/100.0))
except ValueError:
    sys.exit(1)

alltr = dict()
pepgrp = defaultdict(list)
goodtr = set()
ntgs = 0; ngoodtgs = 0;
for row in fdr.get_ids(opts.results):
    if int(row['decoy']) == 1:
	continue
    apex = map(float,row['aggr_Peak_Apex'].split(';'))
    area = map(float,row['aggr_Peak_Area'].split(';'))
    trs = row['aggr_Fragment_Annotation'].split(';')
    rt = float(row['RT'])/60.0
    assay_rt = float(row['assay_rt'])/60.0
    tg = row['transition_group_id']
    score = float(row[scorekey])
    qv = fdr.qvalue(score)
    absint = float(row['Intensity'])
    # print tg,score,qv
    for tr,ap,ar in zip(trs,apex,area):
        alltr[tr] = dict(area=ar,apex=ap,tg=tg,exprt=rt,transnrt=assay_rt,score=score,fdr=qv,intensity=absint)
        if score >= scorethreshold:
	    goodtr.add(tr)
    ntgs += 1
    if score >= scorethreshold:
	ngoodtgs += 1

print "Good TGs:",ngoodtgs,"(Score threshold:",scorethreshold,")"

if ngoodtgs == 0:
    sys.exit(1)

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
	if tr in goodtr:
            pepgrp[(pepid,pz)].append(tr)

lccalfile = opts.results.replace('_table.tsv','_cal.trafoXML')
if os.path.exists(lccalfile):
    for l in open(lccalfile):
	if "<Transformation " in l:
	    m = re.search(r' nrtslope="([^"]*)"',l)
	    nrtslope = float(m.group(1))
	    m = re.search(r' nrtintercept="([^"]*)"',l)
	    nrtintercept = float(m.group(1))
	    break
    print "Extracted from trafoXML file: rt = %s*nrt+%s"%(nrtslope,nrtintercept)
    slope = float(1.0/nrtslope)
    intercept = float(-nrtintercept/nrtslope)
    for tr in goodtr:
	alltr[tr]['expnrt'] = slope*alltr[tr]['exprt']+intercept

elif ntgs >= 2:
    nrtpairs = []
    seen = set()
    for tr in alltr:
	if alltr[tr]['tg'] in seen:
	    continue
	seen.add(alltr[tr]['tg'])
        x = alltr[tr]['transnrt']
        y = alltr[tr]['nrt']
        nrtpairs.append((x,y))
    regress = SLR()
    params = regress.fit(nrtpairs)
    # print "OSWNRT = %f * NRT %+f"%(1/params['slope'],-params['intercept']/params['slope'])
    nrtslope = float(1.0/params['slope'])
    nrtintercept = float(-params['intercept']/params['slope'])
    print "Estimated from results table: rt = %s*nrt+%s"%(nrtslope,nrtintercept)
    for tr in goodtr:
	alltr[tr]['expnrt'] = params['slope']*alltr[tr]['exprt']+params['intercept']

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
        if ele.attrib['id'] not in goodtr:
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

os.makedirs(opts.outdir)

for i,(pepid,charge) in enumerate(pepgrp):
    fn = '%s.%s.50.json'%(pepid,charge)
    data = {}
    data["trlib_z1"] = charge
    data["trlib_pepid"] = pepid
    data["extraction"] = "OpenSWATH"
    data["lccalibration"] = "%s:%s"%(nrtslope,nrtintercept)
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
        if "nrt" not in data and "expnrt" in tr:
            data["nrt"] = tr["expnrt"]
        if "rt" not in data:
            data["rt"] = tr["exprt"]
        if "trlib_mz1" not in data:
            data["trlib_mz1"] = tr["pmz"]
	if "score" not in data:
	    data["score"] = tr["score"]
	if "fdr" not in data:
	    data["fdr"] = tr["fdr"]
	if "intensity" not in data:
	    data["intensity"] = tr["intensity"]
        trdata["name"] = trid
        data["series"].append(trdata)
    # print i,pepid,charge
    wh = open(os.path.join(opts.outdir,fn),'w')
    wh.write(json.dumps(data))
    wh.close()

print "%d peptide groups extracted"%(len(pepgrp),)
