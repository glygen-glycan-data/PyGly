
from . ReferenceTable import ReferenceTable
from . Monosaccharide import constantLookup

class SymbolsTable(ReferenceTable):
    def parseSection(self,name,kv):
        supcls = tuple(sorted(map(constantLookup,kv.get("SuperClass","").split())))
        stem = tuple(sorted(map(constantLookup,kv.get("Stem","").split())))
        mods = []
        for m in kv.get("Mods","").split():
            if ":" in m:
                m = m.split(":")
                m = (int(m[0]),constantLookup(m[1]))
            else:
                m = constantLookup(m)
            mods.append(m)
        mods = tuple(sorted(mods))
        subst = []
        for sub in kv.get("Subst","").split():
            if ":" in sub:
                sub = sub.split(":")
                sub = (int(sub[0]),constantLookup(sub[1]))
            else:
                sub = constantLookup(sub)
            subst.append(sub)
        subst = tuple(sorted(subst))
        key = (supcls,stem,mods,subst)
        return [(key,kv)]
