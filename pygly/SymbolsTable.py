
from ReferenceTable import ReferenceTable
from Monosaccharide import constantLookup

class SymbolsTable(ReferenceTable):
    def parseSection(self,name,kv):
        supcls = tuple(sorted(map(constantLookup,kv["SuperClass"].split())))
        stem = tuple(sorted(map(constantLookup,kv.get("Stem","").split())))
        mods = tuple(sorted(map(constantLookup,kv.get("Mods","").split())))
        subst = tuple(sorted(map(constantLookup,kv.get("Subst","").split())))
        key = (supcls,stem,mods,subst)
        return [(key,kv)]
