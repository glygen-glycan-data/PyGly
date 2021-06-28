#!/bin/env python2

from getwiki import GlycanData
import sys

w = GlycanData()
w.loadsite(sys.argv[1])
