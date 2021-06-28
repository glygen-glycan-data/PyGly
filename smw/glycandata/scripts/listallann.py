#!/bin/env python2

import re, sys
from getwiki import GlycanDataWikiNew

w = GlycanDataWikiNew()
for an in w.iterpages(regex=r'^G.......\.'):
    print an.name
