#!/bin/env python27

from getwiki import GPTWiki
import re

w = GPTWiki()

monos = "NHSF"
for gly in w.iterglycans():
    gsym = gly.get('sym')
    mcnt = {}
    for mono in monos:                                                                                                 
        mcnt[mono] = 0                                                                                                 
        m = re.search(mono+r'(\d+)',gsym)                                                                              
        if m:                                                                                                          
            mcnt[mono] = int(m.group(1))                                                                               
        elif mono in gsym:                                                                                             
            mcnt[mono] = 1                                                                                             
    gly.set('nneuac',mcnt['S'])
    if w.put(gly):
	print gly.get('id')

for pep in w.iterpeptides():
    pepname = pep.get('name')
    pep.set('nox',pepname.count('[Ox]'))
    if w.put(pep):
	print pep.get('id')
