
from ReferenceTable import ReferenceTable
from Monosaccharide import constantLookup

class SymbolsTable(ReferenceTable):
    def parseSection(self,name,kv):
        stem = tuple(sorted(map(constantLookup,kv["Stem"].split())))
        mods = tuple(sorted(map(constantLookup,kv.get("Mods","").split())))
        subst = tuple(sorted(map(constantLookup,kv.get("Subst","").split())))
        key = (stem,mods,subst)
        return [(key,kv)]
