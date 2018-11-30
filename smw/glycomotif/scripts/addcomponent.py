#!/bin/env python27

import os
import sys

pp = os.path.dirname(os.path.abspath(__file__))
component_path = os.path.join(pp, sys.argv[1], "mediawiki", "components.js.txt")

nonredonly = open(os.path.join(pp, "nonredonly.json")).read().strip()
redonly = open(os.path.join(pp, "redonly.json")).read().strip()
topology = open(os.path.join(pp, "topology.json")).read().strip()

component_content = "var topologyComponents = %s\n\n\n\n\n\n\n\nvar redonlyComponents = %s\n\n\n\n\n\n\n\nvar nonredonlyComponents = %s\n\n\n\n\n\n\n\n" %(topology, redonly, nonredonly)

f = open(component_path, "w")
f.write(component_content)
f.close()
