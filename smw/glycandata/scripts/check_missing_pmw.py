#!/bin/env python27

import sys

from getwiki import GlycanDataWiki
w = GlycanDataWiki()

import findpygly
from pygly.CompositionTable import PermethylCompositionTable
from pygly.GlycanFormatter import GlycoCTFormat

pctable = PermethylCompositionTable()
glycoctformat = GlycoCTFormat()

for g in w.iterglycan():
    if not g.has_annotations(property='PermethylatedMW',source='EdwardsLab'):
	glycan = g.getGlycan()
        if not glycan:
            continue
        for m in glycan.all_nodes():
            try:
                eltcomp = m.composition(pctable)
            except KeyError:
                print g.get('accession'),glycoctformat.mtoStr(m)

