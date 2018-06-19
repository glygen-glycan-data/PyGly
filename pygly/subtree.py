from combinatorics import select
import sys, traceback
from MonoFormatter import LinCodeSym, IUPACSym
from collections import defaultdict
from CompositionTable import ResidueCompositionTable
from ElementMass import MonoisotopicElementMass




class subtreeCompare:

    def __init__(self):
        pass

    def compare(self,G1,G2,verbose=False,visibleonly=True):
        if not G1.undetermined() and not G2.undetermined():
            if verbose:
                print >>sys.stderr, "Comparing:\n",G1.str(),"\nwith\n",G2.str()
            return self.compare_(G1.root(),G2.root(),visibleonly=visibleonly)
        nG1 = G1.instantiation_count();
        nG2 = G2.instantiation_count();
	# print nG1,nG2
        # if not self.instSizeTest(nG1,nG2):
        #     return False
        cmp = {}
        for i1,G1i in enumerate(G1.instantiations()):
            for i2,G2i in enumerate(G2.instantiations()):
                cmp[(i1,i2)] = self.compare_(G1i.root(),G2i.root(),visibleonly=visibleonly)
		if self.allmatch() and not cmp[(i1,i2)]:
		    return False
	if self.allmatch():
	    return True
        for sel in select(range(nG2),nG1):
            all = True
            for i1,i2 in zip(range(nG1),sel):
                if not cmp[(i1,i2)]:
                    all = False
                    break
            if all:
                return True
        return False

    def compare_(self,r1,r2,root=True,visibleonly=False):

        if not self.nodeEqual(r1,r2,root,visibleonly):
            return False

        l1 = r1.links()
        l2 = r2.links()

        if not self.linksSizeTest(len(l1),len(l2)):
            return False

        for sel in select(l2,len(l1)):
	    any = False
            for a,b in zip(l1,sel):
                if not self.linkEqual(a,b):
		    any = True
                    break
	    if any:
		continue
            for a,b in zip(l1,sel):
                if not self.compare_(a.child(),b.child(),root=False,visibleonly=visibleonly):
                    any = True
		    break
	    if any:
		continue
            return True
        return False

class subtreeContainedIn(subtreeCompare):
    def linksSizeTest(self,len1,len2):
        return (len1 <= len2)
    def allmatch(self):
	return True

class stringentSubtreeContainedIn(subtreeContainedIn):

    def nodeEqual(self,n1,n2,root,vo):
        return n1.equals(n2)

    def linkEqual(self,l1,l2):
        return l1.equals(l2)

class linearcodeSubtreeContainedIn(subtreeContainedIn):

    def __init__(self):
        self.fmt = LinCodeSym()

    def nodeEqual(self,n1,n2,root,vo):
	# print n1, n2
        try:
            lc1 = self.fmt.toStr(n1)
        except:
            # traceback.print_exc()
	    # print 'FALSE'
            return False
        try:
            lc2 = self.fmt.toStr(n2)
        except:
            # traceback.print_exc()
	    # print 'FALSE'
            return False
	# print lc1, lc2
        if lc1 == lc2:
	    # print 'TRUE'	
            return True
	# print 'FALSE'
        return False

    def linkEqual(self,l1,l2):
        return True


ctable = ResidueCompositionTable()
elmt = MonoisotopicElementMass()
def isobaric_monosaccharide(n1,n2,eltcomp=ctable,eltmass=elmt,tolerance=1e-5):
    massself = n1.composition(eltcomp).mass(eltmass)
    massm = n2.composition(eltcomp).mass(eltmass)
    return (abs(massself - massm) < tolerance)

class isobaricSubtreeContainedIn(subtreeContainedIn):

    def nodeEqual(self,n1,n2,root,vo):
        if isobaric_monosaccharide(n1,n2):
            return True
        return False
        
    def linkEqual(self,l1,l2):
        return True

class subtreeEquals(subtreeCompare):
    def linksSizeTest(self,len1,len2):
        return (len1 == len2)
    def allmatch(self):
	return True

class isobaricSubtreeEquals(subtreeEquals):
    def nodeEqual(self,n1,n2,root,vo):
        if isobaric_monosaccharide(n1,n2):
            return True
        return False

    def linkEqual(self,l1,l2):
        return True

class isobaricSubtreeCompatible(isobaricSubtreeEquals):
    def allmatch(self):
	return False
    
class linearcodeSubtreeEquals(subtreeEquals):

    def __init__(self):
        self.fmt = LinCodeSym()

    def nodeEqual(self,n1,n2,root,vo):
        try:
            lc1 = self.fmt.toStr(n1)
        except:
            return False
        try:
            lc2 = self.fmt.toStr(n2)
        except:
            return False
        if lc1 == lc2:
            return True
        return False

    def linkEqual(self,l1,l2):
        return True

class iupacSubtreeEquals(subtreeEquals):

    def __init__(self):
        self.fmt = IUPACSym()

    def nodeEqual(self,n1,n2,root,vo):
        try:
            lc1 = self.fmt.toStr(n1)
        except:
            return False
        try:
            lc2 = self.fmt.toStr(n2)
        except:
            return False
        if lc1 == lc2:
            return True
        return False

    def linkEqual(self,l1,l2):
        return True

class stringentSubtreeEquals(subtreeEquals):

    def nodeEqual(self,n1,n2,root,vo):
	return n1.equals(n2)

    def linkEqual(self,l1,l2):
        return l1.equals(l2)

class compatibleSubtreeEquals(subtreeEquals):

    def nodeEqual(self,n1,n2,root,vo):
	return n1.compatiblewith(n2,root,vo)

    def linkEqual(self,l1,l2):
        return l1.compatiblewith(l2)

    def allmatch(self):
	return False

