import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')

from getwiki import GPTWiki
from collections import defaultdict
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import sys

w = GPTWiki()

if len(sys.argv) < 2:  
    print 'please enter the spectra file name(s)'
    exit(1)
spectrafiles = sys.argv[1:]

pepnrtpairs = defaultdict(list)

nrtobsgt0 = {}

for tgpage in w.itertransgroups():

    spectra = tgpage.get('spectra')
    if spectra not in spectrafiles:
	continue 
    tgid = tgpage.get('id')
    pepid = tgpage.get('peptide')
    peppage = w.get(pepid)
    pepnrt = peppage.get('nrt')
    nrtobs = peppage.get('nrtobs')
    peakrt = tgpage.get('prt')

    if peakrt != None and pepnrt != None:
	if nrtobs == '0':
	    continue
	else:
            if spectra not in nrtobsgt0:
                nrtobsgt0[spectra] = 1
            else:
                nrtobsgt0[spectra] += 1

	if spectra not in pepnrtpairs:
	    pepnrtpairs[spectra] = [(float(peakrt),float(pepnrt))]
	else:
	    pepnrtpairs[spectra].append((float(peakrt),float(pepnrt)))

# f = open('../data/'+sys.argv[1][:5]+'.coefficients.txt','w')
f = sys.stdout
# f.write('\n')

for s in sorted(pepnrtpairs):

    x = []
    y = []

    for i,(peakrt,pepnrt) in enumerate(pepnrtpairs[s]):
	if pepnrt != None:
	    x.append(pepnrt)
	    y.append(peakrt)
	else:
	    # print 'pep nrt not exit:',s,' ',peakrt
	    pass
	
    pepnrt_x = np.array(x)
    peakrt_y = np.array(y)

    slope, intercept, r_value, p_value, std_err = stats.linregress(pepnrt_x,peakrt_y)

    plt.plot(pepnrt_x, peakrt_y, '.')
    plt.plot(pepnrt_x, intercept + slope*pepnrt_x, 'r')
    plt.legend()
    plt.title(s)
    plt.xlabel('Peptide NRT')
    plt.ylabel('Peak RT(min)')
    plt.savefig(s+'.coefficients.png')
    plt.show()
    plt.close()
    
    f.write('\t'.join(map(str,['filename','nrtobsgt0','slope','intercept','r_value']))+'\n')
    f.write(s+'\t'+str(nrtobsgt0[s])+'\t'+str(round(slope,4))+'\t'+str(round(intercept,2))+'\t'+str(round(r_value,4))+'\n')

    if raw_input('Write to wiki? ').lower()[0] == "y":
        spectrapage = w.get(s)
        #    # if nrtobsgt0[s] == 1 or slope < 0.37 or slope > 0.41:
        # 	spectrapage.set('nrtslope','')
        # 	spectrapage.set('nrtintercept','')
        # 	if w.put(spectrapage):
        ##  	    # print 'spectrapage is updated:',s
        # 	    pass
        #        continue
        spectrapage.set('nrtslope',float(slope))
        spectrapage.set('nrtintercept',float(intercept))
        if w.put(spectrapage):
            # print 'spectra is updated:',s
	    pass

# f.close()
