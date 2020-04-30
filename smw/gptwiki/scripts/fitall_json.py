import sys
import json
import Peak
from getwiki import GPTWiki

from scipy.optimize import curve_fit
from scipy import asarray as ar,exp
from scipy import stats
import numpy as np
import urllib

def gauss(x,a,x0,sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))

def fitall(mmulist=[],exprtlist=[]):    
    for mmu in mmulist:
	xl = []
	yl = []
    	filename = pepid+'.'+str(z1)+'.'+str(mmu)+'.'+'json'
        onlinefile = 'http://edwardslab.bmcb.georgetown.edu/~nedwards/dropbox/pBYmLSkGeq/'+spectra+'/'+filename
        json_file = urllib.urlopen(onlinefile)
        try:
            data = json.load(json_file)
        except:
            print 'file not found'
            continue

	for exprt in exprtlist:
	    xl = []
	    yl = []
	    xlist1 = []
	    ylist1 = []
	    xlist2 = []
	    ylist2 = []
	    xlist3 = []
	    ylist3 = []
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
            print 'window=1.0:'

	    if len(xl) < 8:
                xc, yc = Peak.centroid(xlist2,ylist2)
                xl, yl = Peak.get_startpoint(xc, yc)
                print 'window=2.0:'
	        if len(xl) < 8:
                    xc, yc = Peak.centroid(xlist3,ylist3)
                    xl, yl = Peak.get_startpoint(xc, yc)
                    print 'window=4.0:'

	    if len(xl) >= 8:
	        x = np.array(xl)
        	y = np.array(yl)
        	mean = sum(x*y)/sum(y)
        	sigma = np.sqrt(sum(y*(x-mean)**2)/sum(y))
        	try:
            	    popt,pcov = curve_fit(gauss,x,y,p0=[max(y),mean,sigma])
        	except:
            	    print 'fitting failed'
            	    continue
        	y2 = gauss(x,*popt)
	        print 'before fitting, exprt=',exprt
	        if popt.any():
        	    height = popt[0]
            	    adjrt = popt[1]
            	    print 'after fitting, peakrt=',adjrt
                    std_dev = popt[2]*6
            	    FWHM = 2*np.sqrt(2*np.log(2))*popt[2]
            	    area = FWHM*height
#            	    print 'popt=',popt
        	if std_dev > 3.5:
            	    continue
#Method 1
#z = stats.linregress(x,y)
	        slope, intercept, r_value, p_value, std_err = stats.linregress(y,y2)
#       	print 'slope=',round(slope,4)
#       	print 'intercept=',round(intercept,2)
        	print 'r_value=',round(r_value,4)
        	outfile.write('\t'.join(map(str,[tgid,spectra,pepid,z1,mmu,exprt,adjrt,r_value,std_dev,height,FWHM,area,'\n'])))
	        if tgid not in updated:
		    tg.set('prt',adjrt)
		    updated.add(tgid)
		    if w.put(tg):
			print 'tg is updated:',tgid
	json_file.close()

if len(sys.argv) < 2:  
    print 'please enter the spectra file name(s)'
    exit(1)
spectrafiles = sys.argv[1:]

outfile = open('../data/'+sys.argv[1][:5]+'.fitall.txt','w')

outfile.write('\t'.join(map(str,['TransGroup','Spectra','PeptideID','PrecZ','mmu','ExpRT','AdjRT','R_Value','Std_Dev','Height','FWHM','Area','\n'])))

w = GPTWiki()

for tg in w.itertransgroups():

    tgid = tg.get('id')
    pepid = tg.get('peptide')
    z1 = tg.get('z1')
    exprt = tg.get('rt')
    spectra = tg.get('spectra')
    if spectra not in spectrafiles:
        continue
    scans = tg.get('scans')
    hcds = []
    HCDs = {}
    for i, s in enumerate(scans):
	try:
	    j = s.index('HCD')
	    hcds.append(float(s[j+1]))
	except:
	    continue

    if max(hcds) - float(exprt) > 4:
	HCDs[tg]=hcds

    adjrt = 0.0
    nrt = 0.0
    popt = []
    std_dev = 0.0
    height = 0.0
    FWHM = 0.0
    area = 0.0
    mmu = [5,10,100]
   
    updated = set() 
    if tg in HCDs:
       	print 'fit all scans:',tgid
	fitall(mmu,HCDs[tg])

    else:
	print 'fit one scan:',tgid
	fitall(mmu,[exprt])

outfile.close()

