#!/bin/env python2

import sys
from collections import defaultdict 

from getwiki import GlycanData, Glycan
w = GlycanData()

property = sys.argv[1]
source = sys.argv[2]

def read_aliases(args):
    if len(args) == 0:
        for it in sys.stdin:
            yield it.split()[:2]
    else:
        for fn in args:
            for it in open(fn):
                yield it.split()[:2]

names = defaultdict(set)
for acc,name in read_aliases(sys.argv[3:]):
    names[acc].add(name)

# print names
# for acc in sorted(names):
#     m = w.get(acc)
#     if not m:
# 	continue
    
for m in w.iterglycan():
    acc = m.get('accession')
    m.delete_annotations(property="Name",type="Name",source=source)
    m.delete_annotations(property=property,type="Name",source=source)
    # try:
    #     thenames = set(m.get_annotation_value(property="Name",source=source,type="Name"))
    # except LookupError:
    # 	thenames = set()
    m.set_annotation(value=list(names[acc]),property=property,source=source,type="Name")
    if w.put(m):
        print acc
