from getwiki import GPTWiki

import re

w = GPTWiki()
trs = {}

outfile = open('../data/alltransitions.tsv','w')
outfile.write('\t'.join(map(str,['trid','tgid','pepid','mw','prod.z','prodmz','nrt','precmz','prec.z','glycan','label','relint','sequence','\n'])))

for tgpage in w.itertransgroups():
    
    tgid = tgpage.get('id')
    pepid = tgpage.get('peptide')
    pepage = w.get(pepid)
    try:
        glycan = re.search('\[H(.)*\]',pepage.get('name')).group(0)[1:-1]
    except:
	continue
    mw = pepage.get('mw')
    sequence = pepage.get('sequence')
    precmz = tgpage.get('mz1')
    nrt = tgpage.get('nrt')
    transitions = tgpage.get('transitions')
    
    if nrt == None:
	continue

    for (trid,relint) in transitions:
	if trid not in trs:
	    trs[trid] = {}
	else:
	    continue	
	trs[trid]['pepid'] = pepid
        trs[trid]['tgid'] = tgid
        trs[trid]['precmz'] = precmz
        trs[trid]['nrt'] = nrt
	trs[trid]['relint'] = relint
	trs[trid]['glycan'] = glycan

	trpage = w.get(trid)
	prodmz = trpage.get('mz2')
	label = trpage.get('label')
	z1 = trpage.get('z1')
	z2 = trpage.get('z2')

	trs[trid]['prodmz'] = prodmz
	trs[trid]['label'] = label
	trs[trid]['mw'] = mw
	trs[trid]['z1'] = z1
	trs[trid]['z2'] = z2
	trs[trid]['sequence'] = sequence
	trs[trid]['label'] = label

	outfile.write('\t'.join(map(str,[trid,trs[trid]['tgid'],trs[trid]['pepid'],trs[trid]['mw'],trs[trid]['z2'],trs[trid]['prodmz'],trs[trid]['nrt'],trs[trid]['precmz'],trs[trid]['z1'],trs[trid]['glycan'],trs[trid]['label'],trs[trid]['relint'],trs[trid]['sequence'],'\n'])))

outfile.close()
	   

