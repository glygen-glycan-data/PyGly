#!/bin/env python2

import re, sys
from getwiki import GlycanData

w = GlycanData()
for g in w.iterglycanid():
    print g
