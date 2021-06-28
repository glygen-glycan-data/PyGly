#!/bin/env python2

import re, sys
from getwiki import GlycanDataWikiNew

w = GlycanDataWikiNew()
for an in w.iterpages(include_categories=['Annotation']):
    print an.name
    w.delete(an.name)
