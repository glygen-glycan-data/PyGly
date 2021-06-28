#!/bin/env python2
import sys
import json
import Peak
from getwiki import GPTWiki

from scipy.optimize import curve_fit
from scipy import asarray as ar,exp
from scipy import stats
import numpy as np
import shutil, os, os.path, csv, re

def gauss(x,a,x0,sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))

def fitall(jsonxicfiles,mmu,initrts):
    allresults = []
    for f,mmu in zip(jsonxicfiles,mmu):
	for rt in initrts:
	    res = fit(f,mmu,rt)
	    if res:
		allresults.append(res)
    if len(allresults) > 0:
	return allresults[0]
    return None

def fit(jsonxicfile,mmu,exprt):
    data = json.loads(open(jsonxicfile).read())
    xl = []; yl = []; xlist1 = []; ylist1 = []; xlist2 = []; ylist2 = []; xlist3 = []; ylist3 = [];
    for p in data['peaks']:
        if (p['rt'] > (float(exprt) - 1.0)) and (p['rt'] < (float(exprt) + 1.0)):
            xlist1.append(p['rt'])
            ylist1.append(p['int'])
        if (p['rt'] > (float(exprt) - 2.0)) and (p['rt'] < (float(exprt) + 2.0)):
            xlist2.append(p['rt'])
            ylist2.append(p['int'])
        if (p['rt'] > (float(exprt) - 4.0)) and (p['rt'] < (float(exprt) + 4.0)):
            xlist3.append(p['rt'])
            ylist3.append(p['int'])

    xc, yc = Peak.centroid(xlist1,ylist1)
    xl, yl = Peak.get_startpoint(xc, yc)
    if len(xl) < 8:
        xc, yc = Peak.centroid(xlist2,ylist2)
        xl, yl = Peak.get_startpoint(xc, yc)
	if len(xl) < 8:
            xc, yc = Peak.centroid(xlist3,ylist3)
            xl, yl = Peak.get_startpoint(xc, yc)
	    if len(xl) < 8:
		return None

    x = np.array(xl)
    y = np.array(yl)
    mean = sum(x*y)/sum(y)
    sigma = np.sqrt(sum(y*(x-mean)**2)/sum(y))
    try:
        popt,pcov = curve_fit(gauss,x,y,p0=[max(y),mean,sigma])
    except:
        return None

    height = float(popt[0])
    adjrt = float(popt[1])
    std_dev = float(popt[2])*6
    FWHM = 2*np.sqrt(2*np.log(2))*float(popt[2])
    area = FWHM*height

    if std_dev > 3.5:
        return None

    y2 = gauss(x,*popt)
    slope, intercept, r_value, p_value, std_err = stats.linregress(y,y2)

    return dict(xicfile=jsonxicfile,mmu=mmu,initrt=exprt,peakrt=adjrt,
                rvalue=r_value,stddev=std_dev,height=height,fwhm=FWHM,
                area=area)

mmus = [5,10,100]
for f in sys.argv[1:]:
    # expect *.SA00000.ME00000.merge.txt files
    xicdir = f.rsplit('.',4)[0]
    base = os.path.split(xicdir)[1]
    peakrt = dict()
    seen = set()
    for row in csv.DictReader(open(f),dialect='excel-tab'):
	gphash = row['GlycopeptideHash']
	charge = int(row['PrecursorCharge'])
	if (gphash,charge) in seen:
	    continue
	seen.add((gphash,charge))
        rt = float(row['RetentionTime'])
        scans = row['Scans']
	rts = [rt]
	scanrts = []
	for sc in scans.split(';'):
	    splscan = re.split(r'[(,]',sc)
	    if splscan[1] == "HCD":
		scanrts.append(float(splscan[2]))
	if abs(min(scanrts)-rt) > 4 or abs(max(scanrts)-rt) > 4:
	    rts = scanrts
        result = fitall(["%s/xic-%s.%d.%d.json"%(xicdir,gphash,charge,mmu) for mmu in mmus],mmus,rts)
	if result:
	    peakrt[(gphash,charge)] = result['peakrt']

    shutil.copyfile(f,f+'.bak')
    writer = None; wh = None
    reader = csv.DictReader(open(f+'.bak'),dialect='excel-tab')
    for row in reader:
	if not writer:
	    headers = reader.fieldnames
	    if "PeakRT" not in headers:
		headers.append("PeakRT")
	    wh = open(f,'w')
	    writer = csv.DictWriter(wh,headers,dialect="excel-tab")
	    writer.writeheader()
	gphash = row['GlycopeptideHash']
	charge = int(row['PrecursorCharge'])
	if (gphash,charge) in peakrt:
	    row['PeakRT'] = peakrt[(gphash,charge)]
	writer.writerow(row)
    if wh:
        wh.close()
