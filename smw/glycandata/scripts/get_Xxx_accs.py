#!/bin/env python27

import sys

from getwiki import GlycanDataWiki
w = GlycanDataWiki()

import findpygly
from pygly.MonoFormatter import IUPACSym
from pygly.GlycanFormatter import GlycoCTFormat

iupacSym = IUPACSym()
glycoctformat = GlycoCTFormat()

for g in w.iterglycan():
    if g.has_annotations(property='XxxCount'):
        xxx_count = int(g.get_annotation_value('XxxCount'))
	if xxx_count > 0:
            glycan = g.getGlycan()
            if not glycan:
                continue
            for m in glycan.all_nodes():
                try:
                    sym = iupacSym.toStr(m)
                    if sym not in ('Man','Gal','Glc','Xyl','Fuc','GlcNAc','GalNAc','NeuAc','NeuGc'):
                        print g.get('accession'),glycoctformat.mtoStr(m)
                except KeyError:
                    print g.get('accession'),glycoctformat.mtoStr(m)
                

        
