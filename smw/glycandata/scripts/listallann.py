#!/bin/env python3.12

import re, sys
from getwiki import GlycanDataWikiNew

w = GlycanDataWikiNew()
for an in w.iterpages(regex=r'^G.......\.'):
    print(an.name)
