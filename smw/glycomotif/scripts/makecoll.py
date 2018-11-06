#!/bin/env python27

from getwiki import GlycoMotifWiki, Collection, GlyTouCanMotif, CCRCMotif, GlycoEpitopeMotif, AllMotif
from getwiki import GlydinMotif, GlydinCummingsMotif, GlydinHayesMotif, GlydinCermavMotif, GlydinSugarbindMotif, GlydinBioligoMotif
import sys

w = GlycoMotifWiki()
current = set()
w.put(Collection(id=GlyTouCanMotif.id,
	         name="GlyTouCan Motifs",
                 contact="Kiyoko Aoki-Kinoshita",
                 email="kkiyoko@soka.ac.jp",
                 url="https://glytoucan.org/Motifs/listAll"))
current.add(GlyTouCanMotif.id)
w.put(Collection(id=CCRCMotif.id,
		 name="CCRC Motifs",
                 contact="Rene Ranzinger",
                 email="rene@ccrc.uga.edu",
                 url="https://github.com/glygen-glycan-data/PyGly/raw/master/smw/glycomotif/data/MotifsMP2.xlsx"))
current.add(CCRCMotif.id)
w.put(Collection(id=GlycoEpitopeMotif.id,
		 name="GlycoEpitope Epitopes",
                 contact="Toshisuke Kawasaki",
                 email="tkawasak@fc.ritsumei.ac.jp",
                 url="https://www.glycoepitope.jp"))
current.add(GlycoEpitopeMotif.id)
w.put(Collection(id=AllMotif.id,
		 name="All Motifs",
                 contact="Nathan Edwards",
                 email="nje5@georgetown.edu",
                 url="",
		 primary=False))
current.add(AllMotif.id)
w.put(Collection(id=GlydinMotif.id,
                 name="Glydin",
                 contact="Julien Mariethoz",
                 email="julien.mariethoz@sib.swiss",
                 url="https://github.com/glygen-glycan-data/PyGly/raw/master/smw/glycomotif/data/epitopes.xlsx",
		 primary=False))
current.add(GlydinMotif.id)
w.put(Collection(id=GlydinCummingsMotif.id,
                 name="Glydin - Cummings",
                 contact="Julien Mariethoz",
                 email="julien.mariethoz@sib.swiss",
                 url="https://github.com/glygen-glycan-data/PyGly/raw/master/smw/glycomotif/data/epitopes.xlsx"))
current.add(GlydinCummingsMotif.id)
w.put(Collection(id=GlydinHayesMotif.id,
                 name="Glydin - Hayes",
                 contact="Julien Mariethoz",
                 email="julien.mariethoz@sib.swiss",
                 url="https://github.com/glygen-glycan-data/PyGly/raw/master/smw/glycomotif/data/epitopes.xlsx"))
current.add(GlydinHayesMotif.id)
w.put(Collection(id=GlydinCermavMotif.id,
                 name="Glydin - Cermav",
                 contact="Julien Mariethoz",
                 email="julien.mariethoz@sib.swiss",
                 url="https://github.com/glygen-glycan-data/PyGly/raw/master/smw/glycomotif/data/epitopes.xlsx"))
current.add(GlydinCermavMotif.id)
w.put(Collection(id=GlydinSugarbindMotif.id,
                 name="Glydin - SugarBind",
                 contact="Julien Mariethoz",
                 email="julien.mariethoz@sib.swiss",
                 url="https://github.com/glygen-glycan-data/PyGly/raw/master/smw/glycomotif/data/epitopes.xlsx"))
current.add(GlydinSugarbindMotif.id)
w.put(Collection(id=GlydinBioligoMotif.id,
                 name="Glydin - BiOligo",
                 contact="Julien Mariethoz",
                 email="julien.mariethoz@sib.swiss",
                 url="https://github.com/glygen-glycan-data/PyGly/raw/master/smw/glycomotif/data/epitopes.xlsx"))
current.add(GlydinBioligoMotif.id)

for c in w.itercollection():
    if not c:
	continue
    if c.get('id') not in current:
        print "Deleting:",c.get('pagename')
        w.delete(c.get('pagename'))
