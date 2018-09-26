#!/bin/env python27

from getwiki import GlycoMotifWiki, Collection, GlyTouCanMotif, CCRCMotif, GlycoEpitopeMotif, AllMotif

import sys

w = GlycoMotifWiki()
w.deletemany('Collection',verbose=True)
w.put(Collection(id=GlyTouCanMotif.id,
	         name="GlyTouCan Motifs",
                 contact="Kiyoko Aoki-Kinoshita",
                 email="kkiyoko@soka.ac.jp",
                 url="https://glytoucan.org/Motifs/listAll"))
w.put(Collection(id=CCRCMotif.id,
		 name="CCRC Motifs",
                 contact="Rene Ranzinger",
                 email="rene@ccrc.uga.edu",
                 url=""))
w.put(Collection(id=GlycoEpitopeMotif.id,
		 name="GlycoEpitope Epitopes",
                 contact="Toshisuke Kawasaki",
                 email="tkawasak@fc.ritsumei.ac.jp",
                 url="https://www.glycoepitope.jp"))
w.put(Collection(id=AllMotif.id,
		 name="All Motifs",
                 contact="Nathan Edwards",
                 email="nje5@georgetown.edu",
                 url=""))
