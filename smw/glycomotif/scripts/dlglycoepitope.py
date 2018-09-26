#!/bin/env python27

import sys,os,os.path

from subprocess import Popen, PIPE
def urlopen(url):
    return Popen("wget -q -O - '%s'"%(url,), shell=True, bufsize=4096, stdout=PIPE).stdout

import re, time
lastpage = 0
for l in urlopen("https://www.glycoepitope.jp/epitopes/epitope_list"):
    m = re.search(r'page=(\d+)',l)
    if m and lastpage < int(m.group(1)):
	lastpage = int(m.group(1))
time.sleep(1)


listpages = ["https://www.glycoepitope.jp/epitopes/epitope_list"]
for p in range(2,lastpage+1):
    listpages.append("https://www.glycoepitope.jp/epitopes/epitope_list?page=%d"%(p,))

from collections import defaultdict
epitopes = defaultdict(dict)
for listurl in listpages:
    for l in urlopen(listurl):
	m = re.search(r'<td width="120px"><a href="/epitopes/(\d+)">(EP\d+)</a></td>',l)
	if m:
	     id = int(m.group(1))
	     acc = m.group(2)
	     epitopes[acc]['id'] = id
    time.sleep(1)

for acc in sorted(epitopes):
    for l in urlopen("https://www.glycoepitope.jp/epitopes/%d"%(epitopes[acc]['id'],)):
	m = re.search(r'<th>Epitope ID</th><td>(EP\d+)</td>',l)
	if m:
	    assert m.group(1) == acc
	    epitopes[acc]['acc'] = m.group(1)
	    continue
	m = re.search(r'<th>Epitope name</th><td>(.*)</td>',l)
	if m:
	    epitopes[acc]['name'] = m.group(1)
	    continue
	m = re.search(r'<th>Sequence</th><td>(.*)</td>',l)
	if m:
	    if m.group(1).split()[0] == m.group(1):
	        epitopes[acc]['sequence'] = m.group(1)
	    continue
	m = re.search(r'<a href="https://glytoucan.org/Structures/Glycans/(.*)">',l)
	if m:
	    epitopes[acc]['glytoucan'] = m.group(1)
	    continue
    time.sleep(1)

print "acc\tid\tglytoucan\tsequence\tname"
for acc in sorted(epitopes):
    print "\t".join(map(str,map(lambda h: epitopes[acc].get(h,""),("acc","id","glytoucan","sequence","name"))))

