#!/bin/env python2
from getwiki import GPTWiki
import sys, re

w = GPTWiki()

# if len(sys.argv) < 2:  
#     print 'please enter the spectra file name regex'
#     exit(1)

for spectrapage in w.iterspec(acqtype="DDA"):
  # if not re.search(sys.argv[1],spectrapage.get('name')):
  #   continue
  print spectrapage.get('name')
  nrtslope = spectrapage.get('nrtslope')
  nrtintercept = spectrapage.get('nrtintercept')
  if not nrtslope or not nrtintercept:
    print "No NRT slope or intercept"
    continue
  for tgpage in w.itertgs(spectra=spectrapage.get('name')):
    tgid = tgpage.get('id')
    peakrt = tgpage.get('prt')
    nrt = 0.0

    if peakrt != None:
	nrt = (peakrt - nrtintercept)/nrtslope
	tgpage.set('nrt',nrt)
	if w.put(tgpage):
	    print tgid
    else:  
        tgpage.set('nrt','')
        if w.put(tgpage):
            print tgid
