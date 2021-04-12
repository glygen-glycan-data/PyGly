#!/bin/env python2

from getwiki import GPTWiki

import time, sys

w = GPTWiki()

for sp in w.iterspec(method=sys.argv[1]):
    sp.set('inst',sys.argv[2])
    sp.set('type',sys.argv[3])
    if w.put(sp):
	print sp.get('id')

