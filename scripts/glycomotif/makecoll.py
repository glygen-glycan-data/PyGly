#!/bin/env python27

import findpygly
from pygly.GlycoMotifWiki import GlycoMotifWiki, Collection

w = GlycoMotifWiki()
w.deletemany('Collection')
w.put(Collection(pagename="GlyTouCan",
                 contact="Kiyoko Aoki-Kinoshita",
                 email="kkiyoko@soka.ac.jp",
                 url="https://glytoucan.org/Motifs/listAll"))
w.put(Collection(pagename="UGA-CCRC",
                 contact="Rene Ranzinger",
                 email="rene@ccrc.uga.edu",
                 url=""))
w.put(Collection(pagename="GlyGen",
                 contact="Nathan Edwards",
                 email="nje5@georgetown.edu",
                 url=""))
