#!/bin/env python3.12

import findpygly
from pygly.Glycan import Glycan

print("\t".join(["IUPAC","Subsumption"]))
for sym in sorted(Glycan.iupac_composition_syms + ["Xxx"],key=lambda s: Glycan.subsumption_relationships.get(s,s)):
    if sym in Glycan.subsumption_relationships:
        print("\t".join([sym,Glycan.subsumption_relationships[sym]]))
    else:
        print("\t".join([sym,"-"]))
for sym in Glycan.subst_composition_syms + ["X"]:
    print("\t".join([sym,"-"]))
