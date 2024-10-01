#!/bin/env python3.12

import re, sys
from getwiki import GlycanData

w = GlycanData()
for g in w.iterglycanid():
    print(g)
