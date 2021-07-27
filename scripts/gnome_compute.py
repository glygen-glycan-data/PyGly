#!/bin/env python2

import sys
import findpygly
from pygly.GNOme import main
if len(sys.argv) == 1:
    sys.argv[1:] = "compute -v -v -v -v".split()
main()

