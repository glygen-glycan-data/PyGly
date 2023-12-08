#!/bin/env python2
import sys

from collections import defaultdict

import findpygly

# use the current raw GNOme subsumption dump
# Back-fill with GlyTouCan if needed.

from pygly.GNOme import SubsumptionGraph

gnome = SubsumptionGraph()
gnome.loaddata(sys.argv[1])
sys.argv.pop(1)

accs = set()
for f in sys.argv[1:]:
    accs.update(open(f).read().split())

newaccs0 = set(accs)

while True:
  newaccs = set()

  for acc in newaccs0:
      newaccs.add(acc)
      newaccs.update(filter(lambda anc: anc.startswith('G'),gnome.ancestors(acc)))
      topo = gnome.get_topology(acc)
      comp = gnome.get_composition(acc)
      if comp:
          newaccs.update(gnome.has_composition(comp))
      elif topo:
          newaccs.update(gnome.has_topology(topo))

  if newaccs == newaccs0:
      break

  newaccs0 = set(newaccs)

newaccs = (newaccs-accs)

for acc in sorted(newaccs):
    print acc

