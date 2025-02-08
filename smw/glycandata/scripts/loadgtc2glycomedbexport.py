#!/bin/env python3.12

import sys, re
from collections import defaultdict

from getwiki import GlycanData
w = GlycanData()

import findpygly
from pygly.GlycanResource import GlycomeDBExport

gdbe = GlycomeDBExport()

def validatecfg(s):
    assert re.search(r'^carb([NO]link|Synthe|Glipid)_\d+_(A|P|D000)$',s)
    return s

def validatecarbbank(s):
    assert re.search(r'^\d+[ab]?$',s)
    return s

xrefmap = {'carbbank': ('Carbbank(CCSB)',validatecarbbank),
           'glycosciences.de': ('GLYCOSCIENCES.de',int),
           'cfg': ('CFG',validatecfg),
           'bcsdb': ('BCSDB',int),
           'glycobase(lille)': (None,str),
           'glycomedb': ('GlycomeDB',int)
          }

xref = defaultdict(lambda : defaultdict(set))

for xrefid,gtcacc,resource in gdbe.allgtc():
    try:
        xrefname,valfn = xrefmap[resource]
        xrefid = valfn(xrefid)
    except (LookupError,):
        print("Bad resource:",xrefid,gtcacc,resource,file=sys.stderr)
    except (ValueError,AssertionError):
        print("Bad accession:",xrefid,gtcacc,resource,file=sys.stderr)
    if not xrefname or not xrefid:
        continue
    xref[xrefname][gtcacc].add(xrefid)

for m in w.iterglycan():
    acc = m.get('accession')
    for xrefname in xref:
        if len(xref[xrefname][acc]) > 0:
            m.set_annotation(value=list(xref[xrefname][acc]),property=xrefname,source="GlycomeDB",type="CrossReference")
        else:
            m.delete_annotations(property=xrefname,source="PubChem",type="CrossReference")
    if w.put(m):
        print(acc)
