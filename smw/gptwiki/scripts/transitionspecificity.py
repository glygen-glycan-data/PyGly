import numpy as np
import csv,math
import getwiki
from analysis.fisher import lod, fisher_exact_low, fisher_exact_high

alltrs = {}

prodsd = 0.05
nrtsd = 3
precsd = 30
prodwd = 0.2

def confusekeys(key):
    l = []
    for i in (-1,0,1):
        for j in (-1,0,1):
            for k in (-1,0,1):
                l.append((key[0]+i, key[1]+j, key[2]+k))
    return l

def comparekeys(tr1,tr2):
    key1 = alltrs[tr1]['key']
    key2 = alltrs[tr2]['key']
    confusekeys1 = alltrs[tr1]['confusekeys']
    if key2 in confusekeys1:
        return 1
    else:
	return 0 

def comparenums(tr1,tr2):
    prodmz1 = alltrs[tr1]['prodmz']
    nrt1 = alltrs[tr1]['nrt']
    precmz1 = alltrs[tr1]['precmz']
    
    prodmz2 = alltrs[tr2]['prodmz']
    nrt2 = alltrs[tr2]['nrt']
    precmz2 = alltrs[tr2]['precmz']

    if abs(prodmz2-prodmz1) < prodsd and abs(nrt2-nrt1) < nrtsd and abs(precmz2-precmz1) < precsd:
	return 1
    else:
	return 0

def gettemplist(i,tr):
    l = []
    l.append(tr)
    prodmz1 = alltrs[tr]['prodmz']
    for j in range(i+1,len(alltrs)):
        (tr2,dic2) = sortedlist[j]
	prodmz2 = dic2['prodmz']
	if abs(prodmz2 - prodmz1) < prodwd:
	    l.append(tr2)
	else:
	    break
    for j in range(i-1,0,-1):
        (tr2,dic2) = sortedlist[j]
	prodmz2 = dic2['prodmz']

        if abs(prodmz2 - prodmz1) < prodwd:
            l.append(tr2)
        else:
            break
    return l

def specscore(x,N,n,M):
    if lod(x,N,n,M) < 0:
       # more specific than expected
       sc = -10*math.log(fisher_exact_low(x,N,n,M),10)
    else:
       # less specific than expected
       sc = 10*math.log(fisher_exact_high(x,N,n,M),10)
    return sc

trsfile = csv.DictReader(open('../data/alltransitions.tsv'),delimiter='\t')

for row in trsfile:
    prodmz = float(row['prodmz'])
    nrt = float(row['nrt'])
    precmz = float(row['precmz'])
    trid = row['trid']
    alltrs[trid] = {}
    key = (int(round(float(prodmz)/prodsd)),int(round(float(nrt)/nrtsd)),int(round(float(precmz)/precsd)))
    alltrs[trid]['key'] = key
    alltrs[trid]['confusekeys'] = confusekeys(key)
    alltrs[trid]['glycan'] = row['glycan']
    alltrs[trid]['label'] = row['label']
    alltrs[trid]['prodmz'] = prodmz
    alltrs[trid]['nrt'] = nrt
    alltrs[trid]['precmz'] = precmz
    alltrs[trid]['pepid'] = row['pepid']
    alltrs[trid]['confusecount'] = 0

sortedlist = sorted(alltrs.items(),key=lambda item:item[1]['prodmz'])

templist = []
for i, (tr,dic) in enumerate(sortedlist):
    templist = gettemplist(i,tr)      
    for tr2 in templist:
	if comparekeys(tr,tr2) == 1:
	    if comparenums(tr,tr2) == 1:
	        alltrs[tr]['confusecount'] += 1	    
    templist = []	   
    alltrs[tr]['confusecount'] -= 1      # -1 because it matches to itself

labelspec = {}
nonspeccount = 0
glycanlabelspec = {}

for tr in alltrs:
    glycan = alltrs[tr]['glycan']
    label = alltrs[tr]['label']
    if label[len(label)-1]  == '*':
	label = label[:len(label)-1]
    confusecount = alltrs[tr]['confusecount']
    if label not in labelspec:
	labelspec[label] = {}
	labelspec[label]['non-specificity'] = 0
	labelspec[label]['specificity'] = 0
	if confusecount > 0:
	    labelspec[label]['non-specificity'] += 1
	    nonspeccount += 1
	else:
	    labelspec[label]['specificity'] +=1
    else:
	if confusecount > 0:
            labelspec[label]['non-specificity'] += 1
	    nonspeccount += 1
        else:
            labelspec[label]['specificity'] +=1
    labelspec[label]['total'] = labelspec[label]['non-specificity'] + labelspec[label]['specificity']
    labelspec[label]['percentage'] = labelspec[label]['non-specificity']*100.00/labelspec[label]['total']

    if (glycan,label) not in glycanlabelspec:
        glycanlabelspec[(glycan,label)] = {}
        glycanlabelspec[(glycan,label)]['non-specificity'] = 0
        glycanlabelspec[(glycan,label)]['specificity'] = 0
        if confusecount > 0:
            glycanlabelspec[(glycan,label)]['non-specificity'] += 1
        else:
            glycanlabelspec[(glycan,label)]['specificity'] +=1
    else:
        if confusecount > 0:
            glycanlabelspec[(glycan,label)]['non-specificity'] += 1
        else:
            glycanlabelspec[(glycan,label)]['specificity'] +=1

    glycanlabelspec[(glycan,label)]['total'] = glycanlabelspec[(glycan,label)]['non-specificity'] + glycanlabelspec[(glycan,label)]['specificity']
    glycanlabelspec[(glycan,label)]['percentage'] = glycanlabelspec[(glycan,label)]['non-specificity']*100.00/glycanlabelspec[(glycan,label)]['total']

for label in labelspec:
    if labelspec[label]['total'] >= 10:
        labelspec[label]['specscore'] = specscore(labelspec[label]['non-specificity'],labelspec[label]['total'],nonspeccount,len(alltrs))

sortspecscore = sorted([item for item in labelspec.items() if item[1].has_key('specscore')],
                       key=lambda item:item[1]['specscore'])

labelspecfile = open('../data/labelspecificity.tsv','w')
labelspecfile.write('\t'.join(map(str,['label','nonspec','spec','total','percent','fisherscore']))+'\n')

for (label,dic) in sortspecscore:
    labelspecfile.write('\t'.join(map(str,[label,dic['non-specificity'],dic['specificity'],dic['total'],"%.2f"%dic['percentage'],"%+.2f"%dic['specscore']]))+'\n')

for (glycan,label) in glycanlabelspec:
    if glycanlabelspec[(glycan,label)]['total'] >= 10:
        glycanlabelspec[(glycan,label)]['specscore'] = specscore(glycanlabelspec[(glycan,label)]['non-specificity'],glycanlabelspec[(glycan,label)]['total'],nonspeccount,len(alltrs))

sortglyspecscore = sorted([item for item in glycanlabelspec.items() if item[1].has_key('specscore')],
                          key=lambda item:(item[0][0],item[1]['specscore']))

glylabelspecfile = open('../data/glycanlabelspecificity.tsv','w')
glylabelspecfile.write('\t'.join(map(str,['glycan','label','nonspec','spec','total','percent','fisherscore']))+'\n')

for ((glycan,label),dic) in sortglyspecscore:
    glylabelspecfile.write('\t'.join(map(str,[glycan,label,dic['non-specificity'],dic['specificity'],dic['total'],"%.2f"%dic['percentage'],"%+.2f"%dic['specscore']]))+'\n')

