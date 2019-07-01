#!/bin/env python27

import sys

from getwiki import GlycanData
w = GlycanData()

import findpygly
from pygly.MonoFormatter import IUPACSym
from pygly.GlycanFormatter import GlycoCTFormat

iupacSym = IUPACSym()
glycoctformat = GlycoCTFormat()

for g in w.iterglycan():
    glycan = g.getGlycan()
    if not glycan:
        continue
    for m in glycan.all_nodes():
        try:
            sym = iupacSym.toStr(m)
	    print g.get('accession'),glycoctformat.mtoStr(m),sym
        except KeyError:
            print g.get('accession'),glycoctformat.mtoStr(m),"-"
