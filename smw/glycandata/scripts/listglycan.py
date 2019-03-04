#!/bin/env python27

from getwiki import GlycanDataWiki

w = GlycanDataWiki()

for m in w.iterglycan():
    print m
