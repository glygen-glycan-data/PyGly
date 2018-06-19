
from ElementTable import ElementTable
from IsoShape import IsoShape

class MonoisotopicElementMass(dict):
    def __init__(self):
        elts = ElementTable()
        for e in elts.values():
            self[e.symbol] = e.mono

class AverageElementMass(dict):
    def __init__(self):
        elts = ElementTable()
        for e in elts.values():
            self[e.symbol] = e.ave

class ElementIsotopes(dict):
    def __init__(self):
        elts = ElementTable()
        for e in elts.values():
            self[e.symbol] = e.isotopes
    def cluster(self,comp,maxpos=5):
        sh = IsoShape(self,comp,maxpos)
        return sh.clusterIntensities()

if __name__ == '__main__':
    import sys
    from CompositionTable import Composition
    mt = MonoisotopicElementMass()
    comp = Composition()
    comp.parse(" ".join(sys.argv[1:]))
    print comp.mass(mt)
