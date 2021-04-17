#!/bin/env python2

from getwiki import GPTWiki
import sys

w = GPTWiki()
w.loadsite(sys.argv[1])
