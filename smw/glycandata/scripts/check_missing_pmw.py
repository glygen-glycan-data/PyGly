#!/bin/env python3.12

import sys

from getwiki import GlycanData
w = GlycanData()

import findpygly
from pygly.CompositionTable import PermethylCompositionTable
from pygly.GlycanFormatter import GlycoCTFormat

pctable = PermethylCompositionTable()
glycoctformat = GlycoCTFormat()

for g in w.iterglycan():
    glycan = g.getGlycan()
    if not glycan:
        continue
    for m in glycan.all_nodes():
        try:
            eltcomp = m.composition(pctable)
        except KeyError:
            print(g.get('accession'),glycoctformat.mtoStr(m))
