#!/bin/env python3.12

from getwiki import GlycanData
import sys

w = GlycanData()
w.loadsite(sys.argv[1])
