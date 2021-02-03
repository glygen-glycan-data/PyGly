import json,urllib,re
from getwiki import GPTWiki
import numpy as np
import getwiki
from analysis.fisher import lod, fisher_exact_low, fisher_exact_high
from transitionspecificity import specscore

w = GPTWiki()
trspec = {}
labelspec = {}
glycanlabelspec = {}
smallwindow = 1
largewindow = 6
threshold = 50
nonspeccount = 0
alltrs = set()

for tg in w.itertgs(acqtype='DIA'):
    pepid = tg.get('peptide')
    pepage = w.get(pepid)
    z1 = tg.get('z1')
    spectra = tg.get('spectra')

    if spectra.find('DIA') == -1:
        continue

    try:
        glycan = re.search('\[H(.)*?\]',pepage.get('name')).group(0)[1:-1]
    except:
	continue
    filename = pepid+'.'+str(z1)+'.50.json'
    onlinefile = 'http://edwardslab.bmcb.georgetown.edu/~nedwards/dropbox/pBYmLSkGeq/'+spectra+'/'+filename

    json_file = urllib.urlopen(onlinefile)
    try:
        data = json.load(json_file)
    except:
        print 'file not found'
        continue
    rt = round(data['rt'],3)
    trs = {}

    for tr in data['series']:
        trid = tr['name'].split('_')[2]
	if trid not in alltrs:
	    alltrs.add(trid)
	trpage = w.get(trid)
	label = trpage.get('label')
	if label[len(label)-1] == '*':
	    label = label[:len(label)-1]
        relint = tr['obs_relint']
        pairs = tr['pairs']
        highestinside = 0.0
        highestoutside = 0.0
        for i,(t,h) in enumerate(tr['pairs']):
	    if i == 0:
	        highestinside = h
	        highestoutside = h
	    elif abs(round(t,3) - rt) < smallwindow and h > highestinside:
	        highestinside = h
	    elif abs(round(t,3) - rt) >= smallwindow and abs(round(t,3) - rt) < largewindow and h > highestoutside:
	        highestoutside = h
            else:
	        continue
       
	if highestinside != 0.0:
	    trs[trid] = {}    
	    percent = round(highestoutside*100.0/highestinside,2)
	    trs[trid]['relint'] = relint
	    trs[trid]['percent'] = percent
	    trs[trid]['label'] = label
	    trs[trid]['absint'] = highestinside
	    trs[trid]['glycan'] = glycan

    for tr in trs:    
	label = trs[tr]['label']
	percent = trs[tr]['percent']
	absint = trs[tr]['absint']
	glycan = trs[tr]['glycan']

	if label not in labelspec:
	    labelspec[label] = {}
	    labelspec[label]['absint'] = [absint]
	    labelspec[label]['specificity'] = 0
	    labelspec[label]['non-specificity'] = 0
	    labelspec[label]['percent'] = [percent]
	else:
	    labelspec[label]['percent'].append(percent)	
	    labelspec[label]['absint'].append(absint)
	
	if (glycan,label) not in glycanlabelspec:
            glycanlabelspec[(glycan,label)] = {}
            glycanlabelspec[(glycan,label)]['absint'] = [absint]
            glycanlabelspec[(glycan,label)]['specificity'] = 0
            glycanlabelspec[(glycan,label)]['non-specificity'] = 0
            glycanlabelspec[(glycan,label)]['percent'] = [percent]
        else:
            glycanlabelspec[(glycan,label)]['percent'].append(percent)
            glycanlabelspec[(glycan,label)]['absint'].append(absint)

        if percent < threshold:
	    labelspec[label]['specificity'] += 1
            glycanlabelspec[(glycan,label)]['specificity'] += 1
        else:
	    labelspec[label]['non-specificity'] += 1
            glycanlabelspec[(glycan,label)]['non-specificity'] += 1
	    nonspeccount += 1

for label in labelspec:
    labelspec[label]['total'] = labelspec[label]['specificity'] + labelspec[label]['non-specificity']
    labelspec[label]['percentage'] = labelspec[label]['non-specificity']*100.00/labelspec[label]['total']

    labelspec[label]['absint_median'] = np.median(labelspec[label]['absint'])
    if labelspec[label]['total'] >= 10:
        labelspec[label]['specscore'] = specscore(labelspec[label]['non-specificity'],labelspec[label]['total'],nonspeccount,len(alltrs))

for (glycan,label) in glycanlabelspec:
    glycanlabelspec[(glycan,label)]['total'] = glycanlabelspec[(glycan,label)]['specificity'] + glycanlabelspec[(glycan,label)]['non-specificity']
    glycanlabelspec[(glycan,label)]['percentage'] = glycanlabelspec[(glycan,label)]['non-specificity']*100.00/glycanlabelspec[(glycan,label)]['total']

    glycanlabelspec[(glycan,label)]['absint_median'] = np.median(glycanlabelspec[(glycan,label)]['absint'])
    if glycanlabelspec[(glycan,label)]['total'] >= 10:
        glycanlabelspec[(glycan,label)]['specscore'] = specscore(glycanlabelspec[(glycan,label)]['non-specificity'],glycanlabelspec[(glycan,label)]['total'],nonspeccount,len(alltrs))

sortspecscore = sorted([item for item in labelspec.items() if item[1].has_key('specscore')], key=lambda item:item[1]['specscore'])

sortglycanspecscore = sorted([item for item in glycanlabelspec.items() if item[1].has_key('specscore')],
                          key=lambda item:(item[0][0],item[1]['specscore']))

labelspecfile = open('../data/labelspecificitydia.tsv','w')
labelspecfile.write('\t'.join(map(str,['label','nonspec','spec','total','percent','fisherscore','absint_median']))+'\n')

glycanlabelspecfile = open('../data/glycanlabelspecificitydia.tsv','w')
glycanlabelspecfile.write('\t'.join(map(str,['glycan','label','nonspec','spec','total','percent','fisherscore','absint_median']))+'\n')

for (label,dic) in sortspecscore:
    labelspecfile.write('\t'.join(map(str,[label,dic['non-specificity'],dic['specificity'],dic['total'],"%.2f"%dic['percentage'],"%+.2f"%dic['specscore'],"%+.2f"%dic['absint_median']]))+'\n')

for ((glycan,label),dic) in sortglycanspecscore:
    glycanlabelspecfile.write('\t'.join(map(str,[glycan,label,dic['non-specificity'],dic['specificity'],dic['total'],"%.2f"%dic['percentage'],"%+.2f"%dic['specscore'],"%+.2f"%dic['absint_median']]))+'\n')

