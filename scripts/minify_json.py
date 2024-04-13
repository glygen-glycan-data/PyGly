#!/bin/env python3

import sys,json

data = json.load(open(sys.argv[1]))
json.dump(data,sys.stdout,separators=(',', ':'))
