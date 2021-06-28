
from . ReferenceTable import ReferenceTable
from . Monosaccharide import constantLookup

class ConstantsTable(ReferenceTable):
    def parseSection(self,name,kv):
        const = constantLookup(name)
        return [(const,kv)]
