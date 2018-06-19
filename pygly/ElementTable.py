
from ReferenceTable import ReferenceTable

class Element:
    def __init__(self,name,symbol,mono,ave,isotopes):
        self.name = name
        self.symbol = symbol
        self.mono = mono
        self.ave = ave
        self.isotopes = None
        if isotopes:
            self.isotopes = map(lambda t: tuple(map(float,t.split(':'))),isotopes.split())
            if len(self.isotopes) < 1:
                self.isotopes = None

class ElementTable(ReferenceTable):
    def parseSection(self,name,kv):
        e = Element(name,kv['Symbol'],
                    float(kv['Monoisotopic']),
                    float(kv['Average']),
                    kv.get('Isotopes'))
        return [(name,e)]
