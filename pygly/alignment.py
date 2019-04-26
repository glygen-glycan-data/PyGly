
from Monosaccharide import Linkage, Substituent, Monosaccharide, Mod, Anomer, Stem
from combinatorics import itermatchings, iterecmatchings, itergenmatchings

import inspect
import sys, os.path

def lineno(msg=None):
  callerframerecord = inspect.stack()[1]    # 0 represents this line
                                            # 1 represents line at caller
  frame = callerframerecord[0]
  info = inspect.getframeinfo(frame)
  return 
  if msg:
      print "[%s] %s:%s: %s"%(info.function,os.path.split(info.filename)[1],info.lineno,msg)
  else:
      print "[%s] %s:%s"%(info.function,os.path.split(info.filename)[1],info.lineno)

class Comparitor(object):

    def __init__(self,**kw):
        self._debug = kw.get('debug',False)
        self._verbose = kw.get('verbose',False)

    def debug(self):
        return self._debug

    def verbose(self):
        return self._verbose

    # returns True if a should be considered equivalent to b
    # should be symmetric, reflexive, and transitive

    def eq(self,a,b):
        raise NotImplemented()

    # returns True if a should be considered less than or equal to b
    # in the subset/subsumption/etc. sense.  should be like a partial
    # order, but permit equivalent items. If leq(a,b) and leq(b,a),
    # then eq(a,b). Reflexive, transitive. if not eq(a,b), then
    # anti-symmetric.
    
    def leq(self,a,b):
        raise NotImplemented()

class MonosaccharideComparitor(Comparitor):
    def __init__(self,substcmp=None,sublinkcmp=None,**kw):
        self._substcmp = substcmp
        self._sublinkcmp = sublinkcmp
        super(MonosaccharideComparitor,self).__init__(**kw)

    def substeq(self,a,b):
        return self._substcmp.eq(a,b)

    def sublinkeq(self,a,b):
        return self._sublinkcmp.eq(a,b)

    def substleq(self,a,b):
        return self._substcmp.leq(a,b)

    def sublinkleq(self,a,b):
        return self._sublinkcmp.leq(a,b)

class SubstituentComparitor(Comparitor):
    pass

class LinkageComparitor(Comparitor):
    pass

class SubLinkageComparitor(LinkageComparitor):
    pass

class GlycanEquivalence(Comparitor):

    ### Assumes glycans have the same topology...

    def __init__(self,monocmp=None,linkcmp=None,rootmonocmp=None,**kw):
        self._monocmp = monocmp
	if rootmonocmp:
	    self._rootmonocmp = rootmonocmp
	else:
	    self._rootmonocmp = monocmp
        self._linkcmp = linkcmp
        super(GlycanEquivalence,self).__init__(**kw)

    def rootmonoeq(self,a,b):
        return self._rootmonocmp.eq(a,b)

    def monoeq(self,a,b):
        return self._monocmp.eq(a,b)

    def linkeq(self,a,b):
        return self._linkcmp.eq(a,b)

    def subtree_eq(self,a,b,root=True,mapids=False):

        if root:
            if not self.rootmonoeq(a,b):
                return False
        else:
            if not self.monoeq(a,b):
                return False

        if mapids:
            b.set_id(a.id())

        for ii,jj in itermatchings(a.links(),b.links(),
                                   lambda i,j: self.linkeq(i,j) and self.subtree_eq(i.child(),j.child(),root=False,mapids=mapids)):
            return True

        if mapids:
            b.unset_id()
        
        return False

    def monosaccharide_match(self,a,b):

        if not self.monoeq(a,b):
            return False
	if len(a.parent_links()) != len(b.parent_links()):
	    return False
	if len(a.links(instantiated_only=True)) != len(b.links(instantiated_only=True)):
	    return False
	if len(a.links(instantiated_only=False)) != len(b.links(instantiated_only=False)):
	    return False
        child_links_match = False
        for ii,jj in itermatchings(a.links(instantiated_only=False),b.links(instantiated_only=False),
                                   lambda i,j: self.linkeq(i,j) and self.monoeq(i.child(),j.child())):
            child_links_match = True
            break
        return child_links_match

    def eq(self,a,b):
      
        a.set_ids()        
        b.unset_ids()

        if a.has_root() and b.has_root():
            if not a.undetermined() and not b.undetermined():
                # Simple topologically determined glycan
                return self.subtree_eq(a.root(),b.root(),mapids=True)
            if not self.subtree_eq(a.root(),b.root(),mapids=False):
                # non-composition, but might be undetermined toplogy. Determined part should match.
                return False

        b.set_ids()

        nodeset1 = list(a.all_nodes(subst=False))
        nodeset2 = list(b.all_nodes(subst=False))

        if len(nodeset1) != len(nodeset2):
            return False

        if a.has_root() != b.has_root():
            return False

        # if there are roots, their equivalence has been tested in subtree_align using self.rootmonoeq
        if a.has_root():
            nodeset1.remove(a.root())
        if b.has_root():
            nodeset2.remove(b.root())

        # we assume each node pairs has a single link between it
        linkset1 = dict()
        for l in a.all_links(uninstantiated=True):
            key = (l.parent().id(),l.child().id())
            assert (key not in linkset1)
            linkset1[key] = l

        linkset2 = dict()
        for l in b.all_links(uninstantiated=True):
            key = (l.parent().id(),l.child().id())
            assert (key not in linkset2)
            linkset2[key] = l

        if len(linkset1) != len(linkset2):
            return False
        
        iters = 0
        for ii,jj in iterecmatchings(nodeset1, nodeset2, self.monosaccharide_match):

            iters += 1
            matching = dict(zip(map(lambda m: m.id(),ii),map(lambda m: m.id(),jj)))
            if a.has_root() and b.has_root():
                matching[a.root().id()] = b.root().id()
            # print >>sys.stderr, " ".join(map(lambda i: "%2s"%(i.id(),),ii))
            # print >>sys.stderr, " ".join(map(lambda i: "%2s"%(i.id(),),jj))

            good = True
            for (f,t),l1 in linkset1.items():
                l2 = linkset2.get((matching[f],matching[t]))
                if not l2 or not self.linkeq(l1,l2):
                    good = False
                    break
            if good:
                # print >>sys.stderr, "%d iterations to find an isomorphism"%(iters,)
                return True
                
        return False
        
class CompositionEquivalence(Comparitor):

    ### Assumes we only compare monosaccharides...

    def __init__(self,monocmp=None,**kw):
        self._monocmp = monocmp
        super(CompositionEquivalence,self).__init__(**kw)

    def monoeq(self,a,b):
        return self._monocmp.eq(a,b)

    def monosaccharide_match(self,a,b):

        if not self.monoeq(a,b):
            return False
        return True

    def eq(self,a,b):
      
        a.set_ids()        
        b.set_ids()

        nodeset1 = list(a.all_nodes(subst=False))
        nodeset2 = list(b.all_nodes(subst=False))

        if len(nodeset1) != len(nodeset2):
            return False

        for ii,jj in iterecmatchings(nodeset1, nodeset2, self.monosaccharide_match):
            matching = dict(zip(map(lambda m: m.id(),ii),map(lambda m: m.id(),jj)))
            return True

        return False

class GlycanPartialOrder(Comparitor):

    # Only correct for identical topologies, will need to be
    # extended/changed to handle different topologies with the same
    # composition

    def __init__(self,monocmp=None,linkcmp=None,rootmonocmp=None,**kw):
        self._monocmp = monocmp
	if rootmonocmp:
	    self._rootmonocmp = rootmonocmp
	else:
	    self._rootmonocmp = monocmp
        self._linkcmp = linkcmp

        self._glycompcmp = CompositionPartialOrder(monocmp=self._monocmp)
        self._topocmp = GlycanTopoEqual()

        super(GlycanPartialOrder,self).__init__(**kw)

    def compleq(self,a,b):
        return self._glycompcmp.leq(a,b)

    def topoeq(self,a,b):
        return self._topocmp.eq(a,b)

    def rootmonoleq(self,a,b):
        return self._rootmonocmp.leq(a,b)

    def monoleq(self,a,b):
        return self._monocmp.leq(a,b)

    def linkleq(self,a,b):
        return self._linkcmp.leq(a,b)

    def subtree_leq(self,a,b,root=True):

        if root:
            if not self.rootmonoleq(a,b):
                return False
        else:
            if not self.monoleq(a,b):
                return False

        for ii,jj in itermatchings(a.links(),b.links(),
                                   lambda i,j: self.linkleq(i,j) and self.subtree_leq(i.child(),j.child(),root=False)):
            return True

        return False

    def leq(self,a,b):

        lineno("Start of Subsumption leq")

        # a's composition must be subsumed by b's composition
        if not self.compleq(a,b):
            return False

        lineno("Composition is subsumed")

        # if b is a composition, it doesn't matter what a's topology is
        if not b.has_root():
            return True

        lineno("b is not a composition")

        # b is not a composition, so a must also not be a composition
        if not a.has_root():
            return False

        # Both a and b are now rooted topologies, though not nec. topologically determined
        
        lineno()

        if not a.undetermined() and not b.undetermined():
            # Simple topologically determined glycan
            return self.subtree_leq(a.root(),b.root())

        # at this point at both are rooted and least one is undetermined

        a.set_ids()
        b.set_ids()

        lineno()

        nodeset1 = list(a.all_nodes(subst=False))
        nodeset2 = list(b.all_nodes(subst=False))

        if not self.rootmonoleq(a.root(),b.root()):
            return False

        lineno()

        nodeset1.remove(a.root())
        nodeset2.remove(b.root())

        # we assume each node pairs has a single link between it
        linkset1 = dict()
        for l in a.all_links(uninstantiated=True):
            key = (l.parent().id(),l.child().id())
            assert (key not in linkset1)
            linkset1[key] = l

        linkset2 = dict()
        for l in b.all_links(uninstantiated=True):
            key = (l.parent().id(),l.child().id())
            assert (key not in linkset2)
            linkset2[key] = l

        lineno()

        ninst1 = sum(1 for _ in a.all_links())
        ninst2 = sum(1 for _ in b.all_links())
        nuninst1 = len(linkset1)-ninst1
        nuninst2 = len(linkset2)-ninst2

        if ninst1 < ninst2:
            return False
        if nuninst1 > nuninst2:
            return False
        if len(linkset1) > len(linkset2):
            return False

        lineno()

        iters = 0
        for ii,jj in itergenmatchings(nodeset1, nodeset2, self.monoleq):

            iters += 1
            matching = dict(zip(map(lambda m: m.id(),ii),map(lambda m: m.id(),jj)))
            matching[a.root().id()] = b.root().id()

            # all linkset1 links must be matched, and the only
            # unmatched linkset2 links must be undetermined links

            matchedlinkset2links = set()
            good = True
            for (f,t),l1 in linkset1.items():
                l2 = linkset2.get((matching[f],matching[t]))
                if not l2 or not self.linkleq(l1,l2):
                    good = False
                    break
                else:
                    matchedlinkset2links.add((matching[f],matching[t]))

            if good:
                for (f,t),l2 in linkset2.items():
                    if (f,t) not in matchedlinkset2links and not linkset2[(f,t)].undetermined():
                        good = False
                        break
                
            if good:
                # print >>sys.stderr, "%d iterations to find an isomorphism"%(iters,)
                return True
                
        lineno()

        return False

class CompositionPartialOrder(Comparitor):

    ### Assumes we only compare monosaccharides...

    def __init__(self,monocmp=None,**kw):
        self._monocmp = monocmp
        super(CompositionPartialOrder,self).__init__(**kw)

    def monoleq(self,a,b):
        return self._monocmp.leq(a,b)

    def leq(self,a,b):
      
        a.set_ids()        
        b.set_ids()

        lineno()

        nodeset1 = list(a.all_nodes(subst=False))
        nodeset2 = list(b.all_nodes(subst=False))

        lineno()

        if len(nodeset1) != len(nodeset2):
            return False

        lineno()
             
        for ii,jj in itergenmatchings(nodeset1, nodeset2, self.monoleq):
            return True

        lineno()

        return False

class MonosaccharideEqual(MonosaccharideComparitor):

    def eq(self,a,b):
        if a._anomer != b._anomer:
            return False
        if a._config != b._config:
            return False
        if a._stem != b._stem:
            return False
        if a._superclass != b._superclass:
            return False
        if a._ring_start != b._ring_start:
            return False
        if a._ring_end != b._ring_end:
            return False
        if a._mods != b._mods:
            return False
        any = False
        for ii,jj in itermatchings(a.substituent_links(),b.substituent_links(),
                                   lambda i,j: self.sublinkeq(i,j) and self.substeq(i.child(),j.child())):
            any = True
            break
        if not any:
            return False
        return True
      
class MonosaccharideTopoEqual(MonosaccharideComparitor):

    def eq(self,a,b):
        # Remove anomer, anything else?
        if a._config != b._config:
            return False
        if a._stem != b._stem:
            return False
        if a._superclass != b._superclass:
            return False
        if a._ring_start != b._ring_start:
            return False
        if a._ring_end != b._ring_end:
            return False
        # Should mods be checked without positions?
        if a._mods != b._mods:
            return False
        any = False
        for ii,jj in itermatchings(a.substituent_links(),b.substituent_links(),
                                   lambda i,j: self.sublinkeq(i,j) and self.substeq(i.child(),j.child())):
            any = True
            break
        if not any:
            return False
        return True

class MonosaccharideSubsumed(MonosaccharideComparitor):

    @staticmethod
    def _leq_(a,b):
        if b == None:
            return True
        if a == None:
            return False
        if a == b:
            return True

    def leq(self,a,b):
        if not self._leq_(a._anomer,b._anomer):
            return False
        if not self._leq_(a._config,b._config):
            return False
        if not self._leq_(a._stem,b._stem):
            return False
        if not self._leq_(a._superclass,b._superclass):
            return False
        if not self._leq_(a._ring_start,b._ring_start):
            return False
        if not self._leq_(a._ring_end,b._ring_end):
            return False
        if a._mods != b._mods:
            return False
        any = False
        for ii,jj in itermatchings(a.substituent_links(),b.substituent_links(),
                                   lambda i,j: self.sublinkleq(i,j) and self.substeq(i.child(),j.child())):
            any = True
            break
        if not any:
            return False
        return True
      
class MonosaccharideTopoSubsumed(MonosaccharideSubsumed):

    def leq(self,a,b):
        if not self._leq_(a._config,b._config):
            return False
        if not self._leq_(a._stem,b._stem):
            return False
        if not self._leq_(a._superclass,b._superclass):
            return False
        if not self._leq_(a._ring_start,b._ring_start):
            return False
        if not self._leq_(a._ring_end,b._ring_end):
            return False
        if a._mods != b._mods:
            return False
        any = False
        for ii,jj in itermatchings(a.substituent_links(),b.substituent_links(),
                                   lambda i,j: self.sublinkleq(i,j) and self.substeq(i.child(),j.child())):
            any = True
            break
        if not any:
            return False
        return True
      
class SubstituentEqual(SubstituentComparitor):
    def eq(self,a,b):
	if a._sub != b._sub:
            return False
        return True
            
class LinkageEqual(LinkageComparitor):
    def eq(self,a,b):
	if a._parent_type != b._parent_type:
	    return False
	if a._child_type != b._child_type:
	    return False
	if a._child_pos != b._child_pos:
	    return False
        if a._parent_pos != b._parent_pos:
	    return False
        if a._undetermined != b._undetermined:
            return False
	return True

class LinkageTopoEqual(LinkageComparitor):
    def eq(self,a,b):
	if a._parent_type != b._parent_type:
	    return False
	if a._child_type != b._child_type:
	    return False
        # do not check parent_pos or child_pos
        if a._undetermined != b._undetermined:
            return False
	return True

class LinkageSubsumed(LinkageComparitor):

    @staticmethod
    def _leq_(a,b):
        if b == None:
            return True
        if a == None:
            return False
        if a <= b:
            return True
  
    def leq(self,a,b):
        if not self._leq_(a._parent_type,b._parent_type):
	    return False
        if not self._leq_(a._child_type,b._child_type):
	    return False
        if not self._leq_(a._parent_pos,b._parent_pos):
	    return False
        if not self._leq_(a._child_pos,b._child_pos):
	    return False
        if a._undetermined and not b._undetermined:
            return False
	return True

class LinkageTopoSubsumed(LinkageSubsumed):

    def leq(self,a,b):
        if not self._leq_(a._parent_type,b._parent_type):
	    return False
        if not self._leq_(a._child_type,b._child_type):
	    return False
        if a._undetermined and not b._undetermined:
            return False
	return True

class GlycanEqual(GlycanEquivalence):
    def __init__(self,**kw):
        kw['substcmp']=SubstituentEqual(**kw)
        kw['linkcmp']=LinkageEqual(**kw)
        kw['sublinkcmp']=kw['linkcmp']
        # monotest needs a subst test and a sublink test
        kw['monocmp']=MonosaccharideEqual(**kw)
        super(GlycanEqual,self).__init__(**kw)

class GlycanTopoEqual(GlycanEquivalence):
    def __init__(self,**kw):
        kw['substcmp']=SubstituentEqual(**kw)
        kw['linkcmp']=LinkageTopoEqual(**kw)
        kw['sublinkcmp']=kw['linkcmp']
        # monotest needs a subst test and a sublink test
        kw['monocmp']=MonosaccharideTopoEqual(**kw)
        super(GlycanTopoEqual,self).__init__(**kw)

class GlycanCompEqual(CompositionEquivalence):
    def __init__(self,**kw):
        kw['substcmp']=SubstituentEqual(**kw)
        kw['sublinkcmp']=LinkageTopoEqual(**kw)
        # monotest needs a subst test and a sublink test
        kw['monocmp']=MonosaccharideTopoEqual(**kw)
        super(GlycanCompEqual,self).__init__(**kw)

class GlycanSubsumption(GlycanPartialOrder):
    def __init__(self,**kw):
        kw['substcmp']=SubstituentEqual(**kw)
        kw['linkcmp']=LinkageSubsumed(**kw)
        kw['sublinkcmp']=kw['linkcmp']
        # monotest needs a subst test and a sublink test
        kw['monocmp']=MonosaccharideSubsumed(**kw)
        super(GlycanSubsumption,self).__init__(**kw)

class GlycanCompositionSubsumption(CompositionPartialOrder):
    def __init__(self,**kw):
        kw['substcmp']=SubstituentEqual(**kw)
        kw['sublinkcmp']=LinkageTopoSubsumed(**kw)
        # monotest needs a subst test and a sublink test
        kw['monocmp']=MonosaccharideTopoSubsumed(**kw)
        super(GlycanCompositionSubsumption,self).__init__(**kw)


def items():
    any = False
    for f in sys.argv[1:]:
        any = True
        yield f.strip()
    if not any:
        for l in sys.stdin:
            yield l.strip()

def verify(acc1,acc2,test,result,expected):

    if result != expected:
        print "g1(%s), g2(%s): %s = %s (%s expected): Failed!"%(acc1,acc2,test,result,expected)
        # sys.exit(1)
    else:
        print "g1(%s), g2(%s): %s = %s (%s expected)"%(acc1,acc2,test,result,expected)

def subshow(acc1,acc2,test,result):

    if result:
        print "%s <= %s"%(acc1,acc2)

def compare(gd1,gd2):

    acc1,g1,tacc1,cacc1,bcacc1,l1 = map(gd1.get,("acc","glycan","topo","comp","bcomp","level"))
    acc2,g2,tacc2,cacc2,bcacc2,l2 = map(gd2.get,("acc","glycan","topo","comp","bcomp","level"))

    # Equality testing...

    verify(acc1,acc2,"g1.equals(g2)",g1.equals(g2),acc1==acc2)
    verify(acc1,acc2,"geq.eq(g1,g2)",geq.eq(g1,g2),acc1==acc2)

    # Topology Equality testing...

    if tacc1 and tacc2:

      tg1 = topology(g1)
      tg2 = topology(g2)

      verify(acc1,acc2,"topology(g1).equals(topology(g2))",tg1.equals(tg2),tacc1==tacc2)
      verify(acc1,acc2,"gtopoeq.eq(g1,g2)",gtopoeq.eq(g1,g2),tacc1==tacc2)

    # Composition Equality testing...

    if cacc1 and cacc2:

      cg1 = composition(g1)
      cg2 = composition(g2)

      verify(acc1,acc2,"composition(g1).equals(composition(g2))",cg1.equals(cg2),cacc1==cacc2)
      verify(acc1,acc2,"gcompeq.eq(g1,g2)",gcompeq.eq(g1,g2),cacc1==cacc2)
    
def subsumed(gd1,gd2):

    acc1,g1,tacc1,cacc1,bcacc1,l1 = map(gd1.get,("acc","glycan","topo","comp","bcomp","level"))
    acc2,g2,tacc2,cacc2,bcacc2,l2 = map(gd2.get,("acc","glycan","topo","comp","bcomp","level"))

    subshow(acc1,acc2,"g1 <= g2",subsumption.leq(g1,g2))

       
if __name__ == "__main__":

    import time
    start = time.time()


    from manipulation import Topology, Composition
    topology = Topology()
    composition = Composition()
  
    from GlyTouCan import GlyTouCan
    gtc = GlyTouCan(usecache=True)

    from GlycanFormatter import GlycoCTFormat, WURCS20Format, GlycanParseError
    glycoct_format = GlycoCTFormat()
    wurcs_format = WURCS20Format()

    geq = GlycanEqual()
    gtopoeq = GlycanTopoEqual()
    gcompeq = GlycanCompEqual()
    subsumption = GlycanSubsumption()

    allacc = []; allgly = dict()
    for mass in items():
        count = 0
        for acc1 in sorted(list(gtc.hasmass(mass,2))):

            # if acc1 not in ('G29636OK','G68646SF'):
            #     continue

            # if acc1 not in ('G94131HY', 'G54691SH'):
            #     continue

            topoacc = gtc.gettopo(acc1)
            compacc = gtc.getcomp(acc1)
            bcompacc = gtc.getbasecomp(acc1)
            if acc1 == bcompacc:
                level = 'BaseComposition'
		levelsort = 0
            elif acc1 == compacc:
                level = 'Composition'
		levelsort = 1
            elif acc1 == topoacc:
                level = 'Topology'
		levelsort = 2
            else:
                level = 'Saccharide'
		levelsort = 3

            gly = gtc.getGlycan(acc1)
            if not gly:
                continue

            # if level == 'Saccharide':
            #     continue

	    if not bcompacc:
		continue

            # if level == 'Saccharide' and not gly.undetermined():
            #     continue

            # for m in gly.all_nodes():
            #     print m
            
            allgly[acc1] = dict(acc=acc1,
                                bcomp=bcompacc,
                                comp=compacc,
                                topo=topoacc,
                                level=level,
				levelsort=levelsort,
				mass=gtc.getmass(acc1),
                                glycan=gly)

	check = True
	while check:
	  check = False
	  for acc1 in list(allgly):
	    if allgly[acc1]['bcomp'] not in allgly:
		del allgly[acc1]
	        check = True
		continue
	    if allgly[acc1]['comp'] and allgly[acc1]['comp'] not in allgly:
		del allgly[acc1]
	        check = True
		continue
	    if allgly[acc1]['topo'] and allgly[acc1]['topo'] not in allgly:
		del allgly[acc1]
	        check = True
		continue
	    if not allgly[acc1]['comp'] and allgly[acc1]['level'] != 'BaseComposition':
		del allgly[acc1]
	        check = True
		continue
	    if not allgly[acc1]['topo'] and allgly[acc1]['level'] not in ('BaseComposition','Composition'):
		del allgly[acc1]
	        check = True
		continue

	print "# NODES - %d glycans in molecular weight cluster for %s"%(len(allgly),mass)
	for acc1,g in sorted(allgly.items(),key=lambda t: (t[1].get('levelsort',10),t[0])):
            print acc1,g.get('mass'),g.get('level'),g.get('topo'),g.get('comp'),g.get('bcomp'),
	    gly = g.get('glycan')
            print ("COMP " if (not gly.has_root()) else "") + \
                  ("UNDET " if (gly.undetermined() and gly.has_root()) else "") + \
                  ("FULL" if gly.fully_determined() else "")
	sys.stdout.flush()
        
    from collections import defaultdict
    outedges = defaultdict(list)
    inedges = defaultdict(list)
    allacc = sorted(allgly)
    for acc1 in allacc:
        gly1 = allgly[acc1];
        for acc2 in allacc:
            gly2 = allgly[acc2]
            if acc1 != acc2:
                if subsumption.leq(gly1['glycan'],gly2['glycan']):
                    outedges[acc2].append(acc1)
                    inedges[acc1].append(acc2)

    roots = set(outedges)-set(inedges)

    def printtree(root,edges,multiparents=set(),indent=0):
        print "%s%s%s"%(" "*indent,root,"*" if root in multiparents else ""),allgly[root]['level']
        for ch in edges[root]:
            printtree(ch,edges,multiparents,indent+2)

    def prune_edges(edges):
        toremove = set()
        for n1 in list(edges):
            for n2 in edges[n1]:
              for n3 in edges[n2]:
                  assert n3 in edges[n1], ", ".join([n1, n2, n3,":"]+edges[n1]+[":"]+edges[n2])
                  if n3 in edges[n1]:
                      toremove.add((n1,n3))
        for n1,n3 in toremove:
            edges[n1].remove(n3)
        return edges

    outedges = prune_edges(outedges)
    inedges = prune_edges(inedges)

    multiparents = set()
    for n1 in inedges:
        if len(inedges[n1]) > 1:
            multiparents.add(n1)
      
    print "# TREE"
    for r in inedges:
        if len(inedges[r]) == 0:
            printtree(r,outedges,multiparents)
    sys.stdout.flush()

    print "# EDGES"
    for n1 in outedges:
	if len(outedges[n1]) > 0:
	    print "%s:"%(n1,)," ".join(outedges[n1])
    sys.stdout.flush()

    print "# DONE - Elapsed time %s sec."%(time.time()-start,)
