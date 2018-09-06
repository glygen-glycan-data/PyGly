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
                 url="https://edwardslab.bmcb.georgetown.edu/glycomotif/images/0/09/MotifsMP2.xlsx"))
w.put(Collection(pagename="GlyGen",
                 contact="Nathan Edwards",
                 email="nje5@georgetown.edu",
                 url="https://edwardslab.bmcb.georgetown.edu/glycomotif/GlyGen"))
w.put(Collection(pagename="GlycoEpitope",
                 contact="Toshisuke Kawasaki",
                 email="tkawasak@fc.ritsumei.ac.jp",
                 url="https://www.glycoepitope.jp")
