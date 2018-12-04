#!/bin/env python27

import os, sys


def mv2data(file):
    try:
        os.rename(file, "../data/"+file)
    except OSError:
        pass


mv2data("redonly.json")
mv2data("nonredonly.json")
mv2data("topology.json")
