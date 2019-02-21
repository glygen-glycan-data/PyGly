#!/bin/env python27

import itertools
import findpygly
from pygly.GlycanFormatter import GlycoCTFormat, WURCS20Format

class GlycanLinkEqual():

    def __init__(self):
        pass
    
    def get(self, l1, l2):
        if l1.parent_type() != l2.parent_type():
            return False
        if l1.parent_pos() != l2.parent_pos():
            return False
        if l1.child_type() != l2.child_type():
            return False
        if l1.child_pos() != l2.child_pos():
            return False
        return True


class GlycanLinkCompatibleOneway():

    def __init__(self):
        pass
    
    def get(self, l1, l2):
        if l1.parent_type() != l2.parent_type():
            return False
        if l1.child_type() != l2.child_type():
            return False
        if l1.child_pos() != l2.child_pos():
            return False
        pp1 = l1.parent_pos()
        pp2 = l2.parent_pos()
        if pp1 == None:
            pp1 = set()
        if pp2 == None:
            pp2 = set()
        if len(pp1) == 0 or (pp2 <= pp1):
            return True
        return False

class GlycanLinkCompatibleEitherway():

    def __init__(self):
        pass
    
    def get(self, l1, l2):
        if l1.parent_type() != l2.parent_type():
            return False
        if l1.child_type() != l2.child_type():
            return False
        if l1.child_pos() != l2.child_pos():
            return False
        pp1 = l1.parent_pos()
        pp2 = l2.parent_pos()
        if pp1 == None:
            pp1 = set()
        if pp2 == None:
            pp2 = set()
        if len(pp1) == 0 or len(pp2) == 0 or pp1.intersection(pp2):
            return True
        return False


class GlycanLinkTopologicalSameAs():

    def __init__(self):
        pass

    def get(self, l1, l2):
        if l1.parent_type() != l2.parent_type():
            return False
        if l1.child_type() != l2.child_type():
            return False
        if l1.child_pos() != l2.child_pos():
            return False
        return True

class SubstituentEqual():
    
    def __init__(self):
        pass

    substituentLinkCheck = GlycanLinkEqual()
    
    def get(self, m1, m2):

        subs1 = m1._substituent_links
        subs2 = m2._substituent_links
        if len(subs1) != len(subs2):
            return False
        for newsubs1 in itertools.permutations(subs1):
            flagThisRoundSubCheck = True
            for sub2, sub1 in zip(subs2, newsubs1):
                flagThisRoundSubCheck = self.substituentLinkCheck.get(sub1, sub2) and sub2.child().equals(sub1.child())
                if not flagThisRoundSubCheck:
                    break
            if flagThisRoundSubCheck:
                break
        return flagThisRoundSubCheck



class SubstituentCompatibleOneway(SubstituentEqual):

    substituentLinkCheck = GlycanLinkCompatibleOneway()


class SubstituentEqualEitherway(SubstituentEqual):

    substituentLinkCheck = GlycanLinkCompatibleEitherway()

class MonosaccharideEqual():
    
    subCheck = SubstituentEqual()

    def __init__(self):
        pass
    
    def get(self, m1, m2):
        if m1._anomer != m2._anomer:
            return False
        if m1._config != m2._config:
            return False
        if m1._stem != m2._stem:
            return False
        if m1._superclass != m2._superclass:
            return False
        if m1._ring_start != m2._ring_start:
            return False
        if m1._ring_end != m2._ring_end:
            return False
        if m1._mods != m2._mods:
            return False
        
        if not self.subCheck.get(m1,m2):
            return False

        return True
    

class MonosaccharideCompatibleOneway(MonosaccharideEqual):
    
    subCheck = SubstituentCompatibleOneway()

    def get(self, m1, m2):
        if (m1._anomer) and m1._anomer != m2._anomer:
            return False
        if m1._config != None and m1._config != m2._config:
            return False
        if m1._stem != m2._stem:
            return False
        if m1._superclass != m2._superclass:
            return False
        if m1._ring_start != None and m1._ring_start != m2._ring_start:
            return False
        if m1._ring_end != None and m1._ring_end != m2._ring_end:
            return False
        if m1._mods != m2._mods:
            return False

        if not self.subCheck.get(m1, m2):
            return False

        return True
    

class MonosaccharideCompatibleEitherway(MonosaccharideEqual):
    
    subCheck = SubstituentEqualEitherway()
    
    def get(self, m1, m2):
        if (m1._anomer) and (m2._anomer) and m1._anomer != m2._anomer:
            return False
        if m1._config != m2._config:
            return False
        if m1._stem != m2._stem:
            return False
        if m1._superclass != m2._superclass:
            return False
        if m1._ring_start != m2._ring_start:
            return False
        if m1._ring_end != m2._ring_end:
            return False
        if m1._mods != m2._mods:
            return False

        if not self.subCheck.get(m1, m2):
            return False

        return True


class MonosaccharideTopologySameAs(MonosaccharideEqual):
    subCheck = SubstituentEqualEitherway()

    def get(self, m1, m2):
        if m1._config != m2._config:
            return False
        if m1._stem != m2._stem:
            return False
        if m1._superclass != m2._superclass:
            return False
        if m1._ring_start != m2._ring_start:
            return False
        if m1._ring_end != m2._ring_end:
            return False
        if m1._mods != m2._mods:
            return False

        if not self.subCheck.get(m1, m2):
            return False

        return True

class  RootMonosaccharideTopologySameAs(MonosaccharideTopologySameAs):

    def get(self, m1, m2):
        if m1._config != m2._config:
            return False
        if m1._stem != m2._stem:
            return False
        if m1._superclass != m2._superclass:
            return False
        # if m1._ring_start != None and m2._ring_start != None and m1._ring_start != m2._ring_start:
        #     return False
        # if m1._ring_end != None and m2._ring_end != None and m1._ring_end != m2._ring_end:
        #     return False
        if m1._mods != m2._mods:
            return False

        if not self.subCheck.get(m1, m2):
            return False

        return True

class GlycanEqual():
    
    linkCheck = GlycanLinkEqual()
    monoCheck = MonosaccharideEqual()
    rootMonoCheck = MonosaccharideEqual()
    
    def __init__(self):
        pass
    
    def get(self,g1,g2):
        r1 = g1.root()
        r2 = g2.root()
        return self.recursiveComparison(r1,r2,root=True)
    
    def recursiveComparison(self, r1, r2,root=False):
        links1 = list(r1.links())
        links2 = list(r2.links())
        if ((root and self.rootMonoCheck.get(r1,r2)) or (not root and self.monoCheck.get(r1, r2))) and (len(links1) == len(links2)):
            if len(links1) == 0:
                return True
            else:
                link1 = links1[:]
                for link2 in itertools.permutations(links2):
                    for l1, l2 in zip(link1, link2):
                        equallinks = self.linkCheck.get(l1, l2)
                        equalmono = self.recursiveComparison(l1.child(), l2.child())
                        flag = equallinks and equalmono
                        if not flag:
                            break
                    if flag:
                        return True
                return False
        else:
            return False


class GlycanCompatibleOneway(GlycanEqual):
    
    linkCheck = GlycanLinkCompatibleOneway()
    monoCheck = MonosaccharideCompatibleOneway()
    rootMonoCheck = MonosaccharideCompatibleOneway()


class GlycanCompatibleEitherway(GlycanEqual):

    linkCheck = GlycanLinkCompatibleEitherway()
    monoCheck = MonosaccharideCompatibleEitherway()
    rootMonoCheck = MonosaccharideCompatibleEitherway()



class GlycanTopologySameAs(GlycanEqual):
    
    # This class is used to check whether 2 glycans are topology equal.
    
    linkCheck = GlycanLinkTopologicalSameAs()
    monoCheck = MonosaccharideTopologySameAs()
    rootMonoCheck = RootMonosaccharideTopologySameAs()



class MotifSearchStrict():
    
    def __init__(self):
        pass

    linkCheck = GlycanLinkEqual()
    monoCheck = MonosaccharideEqual()
    rootMonoCheck = MonosaccharideEqual()
    
    def get(self, m, g, reonly = False):
        motif_root = m.root()
        allNodes = g.all_nodes()
        
        for node in allNodes:
            if self.recursiveComparison(motif_root, node, root=True):
                return True
        return False
    
    def recursiveComparison(self, motif_root, gq_node, root=False):
        links1 = list(motif_root.links())
        links2 = list(gq_node.links())
        if ((root and self.rootMonoCheck.get(motif_root,gq_node)) or (not root and self.monoCheck.get(motif_root, gq_node))) and (len(links1) <= len(links2)):
            if len(links1) == 0:
                return True
            else:
                link1 = links1[:]
                for link2 in itertools.permutations(links2):
                    for l1, l2 in zip(link1, link2):
                        equallinks = self.linkCheck.get(l1, l2)
                        equalmono = self.recursiveComparison(l1.child(), l2.child())
                        flag = equallinks and equalmono
                        if not flag:
                            break
                    if flag:
                        return True
                return False
        else:
            return False

    def reducingEndOnly(self, m, g):
        motif_root = m.root()
        node = g.root()

        return self.recursiveComparison(motif_root, node)
    
    def nonReducingEndOnly(self, m, g):
        motif_root = m.root()
        if g.has_root():
            allOtherNodes = list(g.all_nodes())[1:]
        else:
            allOtherNodes = list(g.all_nodes())
        for node in allOtherNodes:
            if self.recursiveComparison(motif_root, node):
                return True
        return False

class MotifSearchAllowWildCards(MotifSearchStrict):

    def __init__(self):
        pass

    linkCheck = GlycanLinkCompatibleEitherway()
    monoCheck = MonosaccharideCompatibleEitherway()
    rootMonoCheck = MonosaccharideCompatibleEitherway()
    
class MotifSearchLoose(MotifSearchStrict):

    def __init__(self):
        pass

    linkCheck = GlycanLinkCompatibleEitherway()
    monoCheck = MonosaccharideCompatibleEitherway()
    rootMonoCheck = MonosaccharideCompatibleOneway()

class MotifSearchTopologicalSameAs(MotifSearchStrict):

    def __init__(self):
        pass

    linkCheck = GlycanLinkTopologicalSameAs()
    monoCheck = MonosaccharideTopologySameAs()
    rootMonoCheck = RootMonosaccharideTopologySameAs()


if __name__ == "__main__":
    seq1 = """RES
    1b:x-dgal-HEX-x:x
    2b:a-dgal-HEX-1:5
    LIN
    1:1o(4+1)2d"""

    seq2 = """RES
    1b:x-dglc-HEX-1:5
    2s:n-acetyl
    3b:b-dgal-HEX-x:x
    4b:a-dgal-HEX-1:5
    LIN
    1:1d(2+1)2n
    2:1o(4+1)3d
    3:3o(4+1)4d"""

    wurcsp = WURCS20Format()
    glycoctp = GlycoCTFormat()

    g1 = glycoctp.toGlycan(seq1)
    g2 = glycoctp.toGlycan(seq2)

    mstsa = MotifSearchTopologicalSameAs()
    print mstsa.get(g1, g2)

