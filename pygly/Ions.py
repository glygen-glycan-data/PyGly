from Glycan import Glycan
from Monosaccharide import Monosaccharide

from CompositionTable import ResidueCompositionTable
from ElementMass import MonoisotopicElementMass
from MonoFactory import MonoFactory
from MonoFormatter import MassSym

from Monosaccharide import Linkage
from MonoFormatter import IUPACSym, LinCodeSym

iupacSym = IUPACSym()
lcSym = LinCodeSym()

ct = ResidueCompositionTable()
elmt = MonoisotopicElementMass()
mf = MonoFactory()
ms = MassSym()


class Ions:
    def __init__(self):
        self._BIONS = []
        self.bholder = Monosaccharide()

    def _BIons(self,glycan,node=None,PARENT=None):
        if node == None:
            node = glycan.root()
            self.bholder = None#node.clone()
            PARENT = node.clone()
            
        child_links = sorted(node.links(),key=Linkage.parent_pos,reverse=True)
        for c in child_links:
            PARENT.add_child(c.child(), child_type=c.child_type(),child_pos=c.child_pos(),parent_type=c.parent_type(),parent_pos=c.parent_pos())
            if self.bholder != None:
                self._BIONS.append(self.bholder)
            self.bholder = c.child()
            self._BIons(glycan,c.child(),PARENT)
        return 0

    def BIons(self,glycan):
        self._BIons(glycan)
        return set(self._BIONS)

if __name__ == '__main__':

    from MonoFactory import MonoFactory
    from Monosaccharide import Linkage
    mf = MonoFactory()

    gc1 = mf.new("GlcNAc")
    gc2 = mf.new("GlcNAc")

    gc1.add_child(gc2,parent_pos=4,
                  parent_type=Linkage.oxygenPreserved,
                  child_type=Linkage.oxygenLost)

    m1 = mf.new('bdMan')
    m2 = mf.new('adMan')
    m3 = mf.new('adMan')
    
    gc2.add_child(m1,parent_pos=4,
                  parent_type=Linkage.oxygenPreserved,
                  child_type=Linkage.oxygenLost)
    m1.add_child(m2,parent_pos=3,
                 parent_type=Linkage.oxygenPreserved,
                 child_type=Linkage.oxygenLost)
    m1.add_child(m3,parent_pos=6,
                 parent_type=Linkage.oxygenPreserved,
                 child_type=Linkage.oxygenLost)

    ###
    #m2.add_child(m3,parent_pos=6,
    #             parent_type=Linkage.oxygenPreserved,
    #             child_type=Linkage.oxygenLost)

    gc3 = mf.new('GlcNAc')
    m2.add_child(gc3,parent_pos=2,
                 parent_type=Linkage.oxygenPreserved,
                 child_type=Linkage.oxygenLost)
    gc4 = mf.new('GlcNAc')
    m3.add_child(gc4,parent_pos=2,
                 parent_type=Linkage.oxygenPreserved,
                 child_type=Linkage.oxygenLost)
    
    g1 = mf.new('bdGal')
    gc3.add_child(g1,parent_pos=4,
                  parent_type=Linkage.oxygenPreserved,
                  child_type=Linkage.oxygenLost)
    g2 = mf.new('bdGal')
    gc4.add_child(g2,parent_pos=4,
                  parent_type=Linkage.oxygenPreserved,
                  child_type=Linkage.oxygenLost)

    s1 = mf.new('Neu5Ac')
    g1.add_child(s1,parent_pos=6,child_pos=2,
                 parent_type=Linkage.oxygenPreserved,
                 child_type=Linkage.oxygenLost)
    s2 = mf.new('Neu5Ac')
    g2.add_child(s2,parent_pos=3,child_pos=2,
                 parent_type=Linkage.oxygenPreserved,
                 child_type=Linkage.oxygenLost)
    
    g = Glycan(gc1)
    print g

    print '\n***test BIONS\n'
    i = Ions()
    
    for ion in i.BIons(g):
        m = 0.0
        print Glycan(ion)
        for node in Glycan(ion).all_nodes():
            n = node.composition(ct)
            mw = n.mass(elmt)
            m += mw
        print 'MW:',m,'\n'
