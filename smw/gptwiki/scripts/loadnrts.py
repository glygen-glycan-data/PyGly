from getwiki import GPTWiki
import sys

w = GPTWiki()

if len(sys.argv) < 2:  
    print 'please enter the spectra file name(s)'
    exit(1)
spectrafiles = sys.argv[1:]

for tgpage in w.itertransgroups():

    spectra = tgpage.get('spectra')
    if spectra not in spectrafiles:
	continue 
    spectrapage = w.get(spectra)
    nrtslope = spectrapage.get('nrtslope')
    nrtintercept = spectrapage.get('nrtintercept')
    tgid = tgpage.get('id')
    peakrt = tgpage.get('prt')
    nrt = 0.0

    if peakrt != None and nrtslope != None:
	nrt = (peakrt - nrtintercept)/nrtslope
	tgpage.set('nrt',nrt)
	if w.put(tgpage):
	    print 'tgpage is updated:',tgid
    else:  
        tgpage.set('nrt','')
        if w.put(tgpage):
            print 'tg is updated:',tgid








