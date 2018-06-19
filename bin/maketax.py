#!/usr/bin/env python27

import sys
import taxonomy
t = taxonomy.NCBITaxonomyTree(sys.argv[1])
t.build(force=True)
