#!/bin/env python27

import sys

from getwiki import GlycanDataWiki
w = GlycanDataWiki()

import findpygly
from pygly.CompositionTable import ResidueCompositionTable
from pygly.GlycanFormatter import GlycoCTFormat

ctable = ResidueCompositionTable()
glycoctformat = GlycoCTFormat()

for g in w.iterglycan():
    glycan = g.getGlycan()
    if not glycan:
        continue
    for m in glycan.all_nodes():
        try:
            eltcomp = m.composition(ctable)
        except KeyError:
            print g.get('accession'),glycoctformat.mtoStr(m)
