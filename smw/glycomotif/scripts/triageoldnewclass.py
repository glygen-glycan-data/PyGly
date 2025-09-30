#!/bin/env python3.12

import sys, os, csv
from collections import defaultdict

alltypes = dict()
allacc = set(filter(lambda acc: acc.startswith('G'),open(sys.argv[1]).read().split()))
old = defaultdict(lambda: defaultdict(set))
for l in open(sys.argv[2]):
    sl = l.strip().split('\t')
    old[sl[0]][sl[1]].add(sl[2])
    if (sl[1],sl[2]) not in alltypes:
        newid = len(alltypes)
        alltypes[(sl[1],sl[2])] = newid
        alltypes[newid] = (sl[1],sl[2])

new = defaultdict(lambda: defaultdict(set))
for l in open(sys.argv[3]):
    sl = l.strip().split('\t')
    if sl[1] == "-":
        continue
    new[sl[0]][sl[1]].add(sl[2])
    if (sl[1],sl[2]) not in alltypes:
        newid = len(alltypes)
        alltypes[(sl[1],sl[2])] = newid
        alltypes[newid] = (sl[1],sl[2])

sigs2acc = defaultdict(set)
for acc in allacc:
    oldsig = []
    for lvl,tset in old[acc].items():
        for t in tset:
            oldsig.append(alltypes[(lvl,t)])
    oldsig = tuple(sorted(oldsig))
    newsig = []
    for lvl,tset in new[acc].items():
        for t in tset:
            newsig.append(alltypes[(lvl,t)])
    newsig = tuple(sorted(newsig))
    if oldsig != newsig:
        sigs2acc[(oldsig,newsig)].add(acc)

for oldsig,newsig in sorted(sigs2acc,key=lambda k: -len(sigs2acc[k])):
    print("#!Accessions:",len(sigs2acc[(oldsig,newsig)]))
    alli = sorted(set(alltypes[i] for i in (oldsig + newsig)))
    removes = [ alltypes[i] for i in oldsig if i not in newsig ]
    adds = [ alltypes[i] for i in newsig if i not in oldsig ]
    for t in alli:
        if t in removes:
            print("#-",t[0],t[1])
        elif t in adds:
            print("#+",t[0],t[1])
        else:
            print("# ",t[0],t[1])
    for acc in sorted(sigs2acc[(oldsig,newsig)])[:10]:
        print(acc)
