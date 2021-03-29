#!/bin/env python27
from getwiki import GPTWiki
import sys, re

w = GPTWiki()

if len(sys.argv) < 2:  
    print 'please enter the spectra file name regex'
    exit(1)

for spectrapage in w.iterspec():
  if not re.search(sys.argv[1],spectrapage.get('name')):
    continue
  nrtslope = spectrapage.get('nrtslope')
  nrtintercept = spectrapage.get('nrtintercept')
  if not nrtslope or not nrtintercept:
    continue
  for tgpage in w.itertgs(spectra=spectrapage.get('name')):
    tgid = tgpage.get('id')
    peakrt = tgpage.get('prt')
    nrt = 0.0

    if peakrt != None and nrtslope != None:
	nrt = (peakrt - nrtintercept)/nrtslope
	tgpage.set('nrt',nrt)
	if w.put(tgpage):
	    print tgid
    else:  
        tgpage.set('nrt','')
        if w.put(tgpage):
            print tgid
