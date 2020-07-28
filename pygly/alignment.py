
from Monosaccharide import Linkage, Substituent, Monosaccharide, Mod, Anomer, Stem
from combinatorics import itermatchings, iterecmatchings, itergenmatchings, iterplacements, itergenmaximalmatchings, choose

import inspect
import sys, os.path
import time
import copy
from collections import defaultdict

verbose = False
def lineno(msg=None):
  if not verbose:
      return
  callerframerecord = inspect.stack()[1]    # 0 represents this line
                                            # 1 represents line at caller
  frame = callerframerecord[0]
  info = inspect.getframeinfo(frame)
  if msg:
      print "[%s, %s] %s:%s: %s"%(time.asctime(),info.function,os.path.split(info.filename)[1],info.lineno,msg)
  else:
      print "[%s, %s] %s:%s"%(time.asctime(),info.function,os.path.split(info.filename)[1],info.lineno)

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

class LinkageComparitorBase(Comparitor):
    pass

class SubLinkageComparitor(LinkageComparitorBase):
    pass

class LinkageComparitor(LinkageComparitorBase):

    def __init__(self, substcmp=None, sublinkcmp=None, linkcmp=None, **kw):
        self._substcmp = substcmp
        self._sublinkcmp = sublinkcmp
        self._linkcmp = linkcmp
        super(Comparitor, self).__init__(**kw)

    def eq(self, a, b):
        pa = a.parent()
        pb = b.parent()

        if pa.is_monosaccharide() != pb.is_monosaccharide():
            return False
        elif pa.is_monosaccharide():
            return self._linkcmp.eq(a, b)
        else:
            # substituent in link
            # assume substituent has and only has one parent, which should be the monosaccharide it attachs to
            la_upper = a.parent().parent_links()[0]
            lb_upper = b.parent().parent_links()[0]
            return self._sublinkcmp.eq(la_upper, lb_upper) and self._linkcmp.eq(a, b) and self._substcmp.eq(pa, pb)

    def leq(self, a, b):
        pa = a.parent()
        pb = b.parent()

        if pa.is_monosaccharide() != pb.is_monosaccharide():
            return False
        elif pa.is_monosaccharide():
            return self._linkcmp.leq(a, b)
        else:
            # substituent in link
            # assume substituent has and only has one parent, which should be the monosaccharide it attachs to
            la_upper = a.parent().parent_links()[0]
            lb_upper = b.parent().parent_links()[0]
            return self._sublinkcmp.leq(la_upper, lb_upper) and self._linkcmp.leq(a, b) and self._substcmp.leq(pa, pb)


def _mindistsfromroot(r):
    if r:
        return _mindistfromroot(r,instonly=False),_mindistfromroot(r,instonly=True)
    return {},{}

def _mindistfromroot(r,d=0,dist=None,instonly=True):
    if not dist:
	dist = dict()
    if dist.get(r.id(),1e+20) > d:
	dist[r.id()] = d
        for l in r.links(instantiated_only=instonly):
	    ch = l.child()
	    dist.update(_mindistfromroot(ch,d+1,dist,instonly))
    return dist

class GlycanEquivalence(Comparitor):

    ### Assumes glycans have the same topology...

    def __init__(self,monocmp=None,linkcmp=None,rootmonocmp=None,**kw):
        self._monocmp = monocmp
	if rootmonocmp:
	    self._rootmonocmp = rootmonocmp
	else:
	    self._rootmonocmp = monocmp
        self._linkcmp = linkcmp
	self.adist = None
	self.bdist = None
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

	assert self.adist and self.bdist
	if self.adist[0].get(a.id()) != self.bdist[0].get(b.id()):
            return False
	if self.adist[1].get(a.id()) != self.bdist[1].get(b.id()):
            return False

        child_links_match = False
        for ii,jj in itermatchings(a.links(instantiated_only=False),b.links(instantiated_only=False),
                                   lambda i,j: self.linkeq(i,j) and self.monoeq(i.child(),j.child())):
            child_links_match = True
            break
        return child_links_match

    def eq(self,a,b):
	
	self.adist = None
	self.bdist = None
      
        a.set_ids()        
        b.unset_ids()

	lineno()

        if a.has_root() and b.has_root():
            if not a.undetermined() and not b.undetermined():
                # Simple topologically determined glycan
                return self.subtree_eq(a.root(),b.root(),mapids=True)
            if not self.subtree_eq(a.root(),b.root(),mapids=False):
                # non-composition, but might be undetermined toplogy. Determined part should match.
                return False

	lineno()

        b.set_ids()

        nodeset1 = list(a.all_nodes(subst=False))
        nodeset2 = list(b.all_nodes(subst=False))

	lineno()

        if len(nodeset1) != len(nodeset2):
            return False

	lineno()

        # if len(nodeset1) == 1:
	#     if a.has_root() and b.has_root():
	#         return self.rootmonoeq(nodeset1[0],nodeset2[0])
	#     return self.monoeq(nodeset1[0],nodeset2[0])

	lineno()

        if a.has_root() != b.has_root():
            return False

        # if there are roots, their equivalence has been tested in subtree_align using self.rootmonoeq
        if a.has_root() and b.has_root():
            nodeset1.remove(a.root())
            nodeset2.remove(b.root())

	lineno()

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

        if len(linkset1) != len(linkset2):
            return False

	lineno()

	# compute distances 
	self.adist = _mindistsfromroot(a.root())
	self.bdist = _mindistsfromroot(b.root())

        iters = 0
        for ii,jj in itergenmatchings(nodeset1, nodeset2, self.monosaccharide_match):

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
		self.adist = None
		self.bdist = None
                return True

        lineno()
                
	self.adist = None
	self.bdist = None
        return False
        
class CompositionEquivalence(Comparitor):

    ### Assumes we only compare monosaccharides...

    def __init__(self,monocmp=None,substcmp=None,**kw):
        self._monocmp = monocmp
        self._substcmp = substcmp
        super(CompositionEquivalence,self).__init__(**kw)

    def nodeeq(self,a,b):
	if a.is_monosaccharide() and b.is_monosaccharide():
            return self._monocmp.eq(a,b)
	elif not a.is_monosaccharide() and not b.is_monosaccharide():
	    return self._substcmp.eq(a,b)
	return False

    def eq(self,a,b):
      
        a.set_ids()        
        b.set_ids()

        nodeset1 = list(a.all_nodes(subst=False,undet_subst=True))
        nodeset2 = list(b.all_nodes(subst=False,undet_subst=True))

        if len(nodeset1) != len(nodeset2):
            return False

        for ii,jj in iterecmatchings(nodeset1, nodeset2, self.nodeeq):
            # matching = dict(zip(map(lambda m: m.id(),ii),map(lambda m: m.id(),jj)))
            return True

        return False

class GlycanPartialOrder(Comparitor):

    def __init__(self,monocmp=None,linkcmp=None,rootmonocmp=None,substcmp=None,**kw):
        self._monocmp = monocmp
	if rootmonocmp:
	    self._rootmonocmp = rootmonocmp
	else:
	    self._rootmonocmp = monocmp
        self._linkcmp = linkcmp

        self._glycompcmp = CompositionPartialOrder(monocmp=self._monocmp,substcmp=substcmp)
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

    def monosaccharide_leq(self,a,b):

        if not self.monoleq(a,b):
            return False
	if len(filter(lambda l: not l.undetermined(),a.parent_links())) < len(filter(lambda l: not l.undetermined(),b.parent_links())):
	    return False
	if len(a.parent_links()) > len(b.parent_links()):
	    return False
	if len(a.links(instantiated_only=True)) < len(b.links(instantiated_only=True)):
            return False
	if len(a.links(instantiated_only=False)) > len(b.links(instantiated_only=False)):
	    return False
	assert self.adist and self.bdist
	# Using uninstantiated edges, a cannot be closer to the root
	if self.adist[0].get(a.id(),1e+20) < self.bdist[0].get(b.id(),1e+20):
            return False
	# Using instantiated edges only, b cannot be closer to the root
	if self.adist[1].get(a.id(),1e+20) > self.bdist[1].get(b.id(),1e+20):
            return False
        return True

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

	self.adist = None
	self.bdist = None

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

	self.adist = _mindistsfromroot(a.root())
	self.bdist = _mindistsfromroot(b.root())

        lineno()

        nodeset1 = list(a.all_nodes(subst=False))
        nodeset2 = list(b.all_nodes(subst=False))

        if not self.rootmonoleq(a.root(),b.root()):
            return False

        lineno()

        nodeset1.remove(a.root())
        nodeset2.remove(b.root())

        # we assume each node pair has a single link between it
        linkset1 = dict()
        for l in a.all_links(uninstantiated=True):
            if l.parent().is_monosaccharide():
                key = (l.parent().id(),l.child().id())
            else:
                key = (l.parent().parent_links()[0].parent().id(), l.child().id())
            assert (key not in linkset1)
            linkset1[key] = l

        linkset2 = dict()
        for l in b.all_links(uninstantiated=True):
            if l.parent().is_monosaccharide():
                key = (l.parent().id(), l.child().id())
            else:
                key = (l.parent().parent_links()[0].parent().id(), l.child().id())
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
        for ii,jj in itergenmatchings(nodeset1, nodeset2, self.monosaccharide_leq):

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

        # print >>sys.stderr, "%d iterations to prove no isomorphism exists"%(iters,)
        return False



class ConnectedNodesCache:

    def __init__(self):
        self.data = {}

    def put(self, g):
        for m in g.all_nodes():
            self.data[m] = {
                1: [{m}]
            }

    def update_cache(self, m, size):

        if size in self.data[m]:
            return
        i = max(self.data[m].keys())

        while i < size:
            i += 1

            res = []
            for currentSet in self.data[m][i-1]:
                for res0 in self.connectedNodesPlusOneSimple(currentSet):
                    if res0 not in res:
                        res.append(res0)
            self.data[m][i] = res

    def connectedNodesPlusOneSimple(self, currentSet):
        res = []
        for n in currentSet:
            res += [x.child() for x in n.links(instantiated_only=False)]
        children = filter(lambda x:x not in currentSet, res)
        children = set(children)

        res = []
        for n in children:
            res0 = currentSet.copy()
            res0.add(n)
            res.append(res0)
        return res

    def get(self, m, size):
        size = int(size)
        self.update_cache(m, size)
        return self.data[m][size]





class SubstructureSearch(GlycanPartialOrder):

    def __init__(self, **kw):
        self.connected_nodes_pre_computed = kw.get("connected_nodes_pre_computed", True)
        if self.connected_nodes_pre_computed:
            self.nodes_cache = kw.get("connected_nodes_cache", ConnectedNodesCache())
        super(SubstructureSearch, self).__init__(**kw)

    def check_links(self, motif, motif_node_set, tg_node_set):
        motif_nodes = motif_node_set
        tg_nodes = tg_node_set

        discover = [motif.root()]
        while len(discover) > 0:
            this_mono = discover.pop()

            for l in this_mono.links():
                p = this_mono
                c = l.child()
                discover.append(c)

                tgp = tg_nodes[motif_nodes.index(p)]
                # tgc = tg_nodes[motif_nodes.index(c)]

                found_matched_link = False
                for ll in tgp.links(instantiated_only=False):
                    if ll.child() not in tg_nodes:
                        continue

                    if motif_nodes.index(c) == tg_nodes.index(ll.child()) and self._linkcmp.leq(l, ll):
                        found_matched_link = True
                        break

                if not found_matched_link:
                    return False

        return True

    # For recursive algorithm
    def connectedNodesPlusOne(self, currentSet, candidates):
        res = []
        for c in candidates:
            newCurrentSet = currentSet.copy()
            newCandidatesSet = candidates.copy()

            newCurrentSet.add(c)
            newCandidatesSet.remove(c)

            issueFlag = False
            for cl in c.links(instantiated_only=False):
                toAdd = cl.child()
                if toAdd in newCurrentSet and toAdd in newCandidatesSet:
                    issueFlag = True
                    break
                else:
                    newCandidatesSet.add(toAdd)
            if issueFlag:
                continue

            r = (newCurrentSet, newCandidatesSet)
            if r not in res:
                res.append(r)
                yield r

    def connectedNodes(self, currentSet, candidates, size):
        # set, set, int
        if len(currentSet) == size:
            yield currentSet
            raise StopIteration
        elif len(candidates) == 0:
            raise StopIteration

        seen = []
        for newCurrentSet, newCandidatesSet in self.connectedNodesPlusOne(currentSet, candidates):
            for r in self.connectedNodes(newCurrentSet, newCandidatesSet, size):
                if r not in seen:
                    seen.append(r)
                    yield r

    # Recursive algorithm
    def allConnectedNodesByRootRecursive(self, r, size):
        res = []
        for cn in self.connectedNodes({r},
                                      set([x.child() for x in r.links(instantiated_only=False)]),
                                      size):
            res.append(cn)

        return res



    def allConnectedNodesByRoot(self, r, size):
        if self.connected_nodes_pre_computed:
            return self.nodes_cache.get(r, size)
        else:
            return self.allConnectedNodesByRootRecursive(r, size)

    def subtree_leq(self,m, tg,root=True):

        if root:
            if not self.rootmonoleq(m, tg):
                return False
        else:
            if not self.monoleq(m, tg):
                return False

        for tg_linkset in choose(tg.links(), len(m.links())):
            for ii,jj in itermatchings(m.links(), tg_linkset,
                                       lambda i, j: self.linkleq(i, j) and self.subtree_leq(i.child(), j.child(), root=False)):
                return True

        return False


    def leq(self, m, tg, rootOnly=False, anywhereExceptRoot=False):

        assert (rootOnly and anywhereExceptRoot) is not True

        if not tg.has_root():
            return False

        if len(list(m.all_nodes())) > len(list(tg.all_nodes())):
            return False

        tg_root = tg.root()
        tmp = []

        if rootOnly:
            tmp = [tg_root]
        elif anywhereExceptRoot:
            tmp = filter(lambda n: n != tg_root, tg.all_nodes())
        else:
            tmp = list(tg.all_nodes())

        potential_TG_root = []
        for tg_root in tmp:
            if self.rootmonoleq(m.root(), tg_root):
                potential_TG_root.append(tg_root)

        # Use subtree algorithm to save runtime
        if not m.undetermined():
            for n in potential_TG_root:
                if self.subtree_leq(m.root(), n):
                    return True

            # return False when no match based on subtree algorithm and the glycan is not undetermined.
            if not tg.undetermined():
                return False

        if tg.undetermined() and self.connected_nodes_pre_computed:
            # pre compute the connected node set
            self.nodes_cache.put(tg)


        #i = 0
        for tg_root in potential_TG_root:

            for tg_nodes in self.allConnectedNodesByRoot(tg_root, len(list(m.all_nodes()))):
                # print map(lambda x: x.external_descriptor_id(), tg_nodes),map(lambda x: x.external_descriptor_id(), list(m.all_nodes()))
                #i+=1
                #j=0
                for monoset_m, monoset_tg in itergenmatchings(list(m.all_nodes()), tg_nodes, self.monoleq):
                    # matching = dict(zip(map(lambda m: m.external_descriptor_id(), monoset_m), map(lambda m: m.external_descriptor_id(), monoset_tg)))
                    # print matching
                    #j+=1
                    #print i, j
                    if self.check_links(m, monoset_m, monoset_tg):
                        return True

        return False


class CompositionPartialOrder(Comparitor):

    ### Assumes we only compare nodes to each other...
    ### Complication is in dealing with floating substituents

    def __init__(self,monocmp=None,substcmp=None,**kw):
        self._monocmp = monocmp
        self._substcmp = substcmp
        super(CompositionPartialOrder,self).__init__(**kw)

    def nodeleq(self,a,b):
	if a.is_monosaccharide() and b.is_monosaccharide():
            return self._monocmp.leq(a,b)
	elif not a.is_monosaccharide() and not b.is_monosaccharide():
	    return self._substcmp.leq(a,b)
	return False

    def leq(self,a,b):
      
        a.set_ids()        
        b.set_ids()

        lineno()

        nodeset1 = list(a.all_nodes(subst=False))
        nodeset2 = list(b.all_nodes(subst=False))

        nodeset1all = list(a.all_nodes(subst=False,undet_subst=True))
        nodeset2all = list(b.all_nodes(subst=False,undet_subst=True))

        nodeset1uds = filter(lambda n: not n.is_monosaccharide(),nodeset1all)
        nodeset2uds = filter(lambda n: not n.is_monosaccharide(),nodeset2all)

        lineno()

        if len(nodeset1) != len(nodeset2):
            return False

        lineno()

        # print len(nodeset1),len(nodeset1uds),len(nodeset2),len(nodeset2uds)

        # for m in nodeset1all:
        #     print m

        # for m in nodeset2all:
        #     print m

        # sys.stdout.flush()

	# allow for # floating substituent monosaccharides not to match, then deal with these in loop
	fs1 = len(nodeset1uds)
	fs2 = len(nodeset2uds)

	if fs1 > fs2:
	    return False

        # we deal with three cases, for now...
        # 1. No floating substituents at all (fs1 == fs2 == 0)
        # 2. Same number of floating substituents (fs1 == fs2)
        # 3. All or nothing: instances of a specific substituent are
        #    either all floating or all linked. If all floating on one
        #    and all linked on the other, must be floating for b and
        #    linked for a.
	
	if fs2 == fs1:
            for ii,jj in itergenmatchings(nodeset1all,nodeset2all,self.nodeleq):
		return True

	else:

            allsubst = set()
            subst1link = defaultdict(int)
            subst1float= defaultdict(int)
            for n in nodeset1all:
                if n.is_monosaccharide():
                    for s in n.substituents():
                        subst1link[s.name()] += 1
                        allsubst.add(s.name())
                else:
                    subst1float[n.name()] += 1
                    allsubst.add(n.name())

            subst2link = defaultdict(int)
            subst2float= defaultdict(int)
            for n in nodeset2all:
                if n.is_monosaccharide():
                    for s in n.substituents():
                        subst2link[s.name()] += 1
                        allsubst.add(s.name())
                else:
                    subst2float[n.name()] += 1
                    allsubst.add(n.name())

            # print allsubst
            # print subst1link
            # print subst1float
            # print subst2link
            # print subst2float

            subst2move = []
            bad = False
            for s in allsubst:
                if (subst1link[s]+subst1float[s]) != (subst2link[s]+subst2float[s]):
                    return False
                elif subst2float[s] > 0 and \
                         subst2link[s] == 0 and \
                         subst1float[s] == 0 and \
                         subst2float[s] == subst1link[s]:
                    subst2move.append(s)
                elif subst2float[s] > 0 and \
                         subst2link[s] == 0 and \
                         subst1link[s] == 0 and \
                         subst2float[s] == subst1float[s]:
                    # Do nothing, but not a problem
                    pass
		elif subst2float[s] == 0 and \
			 subst1float[s] == 0 and \
			 subst1link[s] == subst2link[s]:
		    # Do nothing, but not a problem
		    pass
                else:
                    bad = True
                    break

            if bad:
                raise RuntimeError("Cannot resolve composition-partial order, floating substituents too complicated")

            # print subst2move
            
            # strip specific floating substituents from nodes in nodeset1
            newnodeset1 = []
            for m in nodeset1:
                m1 = m.clone()
                for sl in m1.substituent_links():
                    if sl.child().name() in subst2move:
                        m1.remove_substituent_link(sl)
                newnodeset1.append(m1)

            # Now the monosaccharides should match....
            for ii,jj in itergenmatchings(newnodeset1,nodeset2,self.nodeleq):
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

class MonosaccharideEqualWithWURCSCheck(MonosaccharideEqual):

    def eq(self,a,b):
        if a.external_descriptor() != b.external_descriptor():
            return False
        return super(MonosaccharideEqualWithWURCSCheck, self).eq(a, b)
      
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

class MonosaccharideImageEqual(MonosaccharideComparitor):

    def eq(self,a,b):
        if a._stem != b._stem:
            return False
        if a._superclass != b._superclass:
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

class MonosaccharideMotifComparison(MonosaccharideComparitor):

    def leq(self, m, g):
        if m._anomer and m._anomer != g._anomer:
            return False
        if m._config and m._config != g._config:
            return False
        if m._stem and m._stem != g._stem:
            return False
        if m._superclass != g._superclass:
            return False
        if m._ring_start and m._ring_start != g._ring_start:
            return False
        if m._ring_end and m._ring_end != g._ring_end:
            return False
        if m._mods != g._mods:
            return False

        any = False
        for ii, jj in itermatchings(m.substituent_links(), g.substituent_links(),
                                    lambda i, j: self.sublinkeq(i, j) and self.substeq(i.child(), j.child())):
            any = True
            break
        if not any:
            return False

        return True



class MonosaccharideMotifComparisonOptionalSubst(MonosaccharideComparitor):

    def leq(self, m, g):
        if m._anomer and m._anomer != g._anomer:
            return False
        if m._config and m._config != g._config:
            return False
        if m._stem and m._stem != g._stem:
            return False
        if m._superclass != g._superclass:
            return False
        if m._ring_start and m._ring_start != g._ring_start:
            return False
        if m._ring_end and m._ring_end != g._ring_end:
            return False

        gmod = copy.deepcopy(g._mods)
        for mod in gmod:
            if mod[1] == Mod.aldi:
                gmod.remove(mod)

        if m._mods != gmod:
            return False

        any = False

        msl = m.substituent_links()
        gsl = g.substituent_links()
        if len(msl) > len(gsl):
            return False

        gsl_mandatory, gsl_optional = [], []
        for sl in gsl:
            if sl.child()._sub in [Substituent.sulfate, Substituent.phosphate]:
                gsl_optional.append(sl)
            else:
                gsl_mandatory.append(sl)

        for x in choose(gsl_optional, len(msl) - len(gsl_mandatory)):

            for ii, jj in itermatchings(msl, x + gsl_mandatory,
                                        lambda i, j: self.sublinkeq(i, j) and self.substeq(i.child(), j.child())):
                any = True
                break
            if any:
                break
        if not any:
            return False

        return True




      
class SubstituentEqual(SubstituentComparitor):
    def eq(self,a,b):
	if a._sub != b._sub:
            return False
        return True
    def leq(self,a,b):
	if a._sub != b._sub:
            return False
        return True


class SubstituentEqualWithWURCSCheck(SubstituentEqual):
    def eq(self,a,b):
        if a.external_descriptor() != b.external_descriptor():
            return False
        return super(SubstituentEqualWithWURCSCheck, self).eq(a, b)


            
class LinkageEqualSimple(LinkageComparitorBase):
    def eq(self, a, b):
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

class LinkageTopoEqualSimple(LinkageComparitorBase):
    def eq(self,a,b):
	if a._parent_type != b._parent_type:
	    return False
	if a._child_type != b._child_type:
	    return False
        # do not check parent_pos or child_pos
        if a._undetermined != b._undetermined:
            return False
	return True

class LinkageImageEqualSimple(LinkageComparitorBase):
    def eq(self,a,b):
        if a._undetermined != b._undetermined:
            return False
	return True

class LinkageSubsumedSimple(LinkageComparitorBase):

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

class LinkageTopoSubsumedSimple(LinkageSubsumedSimple):

    def leq(self,a,b):
        if not self._leq_(a._parent_type,b._parent_type):
	    return False
        if not self._leq_(a._child_type,b._child_type):
	    return False
        if a._undetermined and not b._undetermined:
            return False
	return True

class LinkageComparitorMotifSimple(LinkageComparitorBase):

    def leq(self, m, g):
        if m._parent_type != g._parent_type:
            return False
        if m._child_type != g._child_type:
            return False
        if m._child_pos and g._child_pos and len(set(m._child_pos).intersection(set(g._child_pos))) == 0:
            return False
        if m._parent_pos and g._parent_pos and len(set(m._parent_pos).intersection(set(g._parent_pos))) == 0:
            return False
        return True

class SubstLinkageComparitorMotifSimple(LinkageEqualSimple):

    def leq(self,a,b):
        return self.eq(a,b)


class LinkageEqual(LinkageComparitor):

    def __init__(self, **kw):
        kw['substcmp'] = SubstituentEqual(**kw)
        kw['linkcmp'] = LinkageEqualSimple(**kw)
        kw['sublinkcmp'] = kw['linkcmp']
        super(LinkageEqual, self).__init__(**kw)



class LinkageTopoEqual(LinkageComparitor):

    def __init__(self, **kw):
        kw['substcmp'] = SubstituentEqual(**kw)
        kw['linkcmp'] = LinkageTopoEqualSimple(**kw)
        kw['sublinkcmp'] = kw['linkcmp']
        super(LinkageTopoEqual, self).__init__(**kw)


class LinkageImageEqual(LinkageComparitor):

    def __init__(self, **kw):
        kw['substcmp'] = SubstituentEqual(**kw)
        kw['linkcmp'] = LinkageImageEqualSimple(**kw)
        kw['sublinkcmp'] = kw['linkcmp']
        super(LinkageImageEqual, self).__init__(**kw)


class LinkageSubsumed(LinkageComparitor):

    def __init__(self, **kw):
        kw['substcmp'] = SubstituentEqual(**kw)
        kw['linkcmp'] = LinkageSubsumedSimple(**kw)
        kw['sublinkcmp'] = kw['linkcmp']
        super(LinkageSubsumed, self).__init__(**kw)


class LinkageTopoSubsumed(LinkageComparitor):

    def __init__(self, **kw):
        kw['substcmp'] = SubstituentEqual(**kw)
        kw['linkcmp'] = LinkageTopoSubsumedSimple(**kw)
        kw['sublinkcmp'] = kw['linkcmp']
        super(LinkageComparitor, self).__init__(**kw)

class LinkageMotifComparitor(LinkageComparitor):

    def __init__(self, **kw):
        kw['substcmp'] = SubstituentEqual(**kw)
        kw['linkcmp'] = LinkageComparitorMotifSimple(**kw)
        kw['sublinkcmp'] = SubstLinkageComparitorMotifSimple(**kw)
        super(LinkageMotifComparitor, self).__init__(**kw)

class GlycanEqual(GlycanEquivalence):
    def __init__(self,**kw):
        kw['substcmp']=SubstituentEqual(**kw)
        kw['linkcmp']=LinkageEqual(**kw)
        kw['sublinkcmp']=kw['linkcmp']
        # monotest needs a subst test and a sublink test
        kw['monocmp']=MonosaccharideEqual(**kw)
        super(GlycanEqual,self).__init__(**kw)

class GlycanEqualWithWURCSCheck(GlycanEquivalence):
    def __init__(self,**kw):
        kw['substcmp']=SubstituentEqualWithWURCSCheck(**kw)
        kw['linkcmp']=LinkageEqual(**kw)
        kw['sublinkcmp']=kw['linkcmp']
        # monotest needs a subst test and a sublink test
        kw['monocmp']=MonosaccharideEqualWithWURCSCheck(**kw)
        super(GlycanEqualWithWURCSCheck,self).__init__(**kw)

class GlycanTopoEqual(GlycanEquivalence):
    def __init__(self,**kw):
        kw['substcmp']=SubstituentEqual(**kw)
        kw['linkcmp']=LinkageTopoEqual(**kw)
        kw['sublinkcmp']=kw['linkcmp']
        # monotest needs a subst test and a sublink test
        kw['monocmp']=MonosaccharideTopoEqual(**kw)
        super(GlycanTopoEqual,self).__init__(**kw)

class GlycanImageEqual(GlycanEquivalence):
    def __init__(self,**kw):
        kw['substcmp']=SubstituentEqual(**kw)
        kw['linkcmp']=LinkageImageEqual(**kw)
        kw['sublinkcmp']=kw['linkcmp']
        # monotest needs a subst test and a sublink test
        kw['monocmp']=MonosaccharideImageEqual(**kw)
        super(GlycanImageEqual,self).__init__(**kw)

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

class GlyTouCanMotif(SubstructureSearch):

    def __init__(self, **kw):
        kw["monocmp"] = MonosaccharideMotifComparison(
            substcmp=SubstituentEqual(),
            sublinkcmp=LinkageEqual()
        )
        kw["linkcmp"] = LinkageMotifComparitor()
        super(GlyTouCanMotif, self).__init__(**kw)



class MotifAllowOptionalSub(SubstructureSearch):

    def __init__(self, **kw):
        kw["monocmp"] = MonosaccharideMotifComparisonOptionalSubst(
            substcmp=SubstituentEqual(),
            sublinkcmp=LinkageEqual()
        )
        kw["linkcmp"] = LinkageMotifComparitor()
        super(MotifAllowOptionalSub, self).__init__(**kw)





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
        print "%s %s %s"%(acc1,test,acc2)
    else:
        print "%s %s/%s %s"%(acc1,test[0],test[-1],acc2)

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

    from collections import defaultdict

    from manipulation import Topology, Composition
  
    # from GlyTouCan import GlyTouCan
    from GlycanResource import GlyTouCan
    gtc = GlyTouCan(usecache=False)

    from GlycanFormatter import GlycoCTFormat, WURCS20Format, GlycanParseError
    glycoct_format = GlycoCTFormat()
    wurcs_format = WURCS20Format()

    geq = GlycanEqual()
    gtopoeq = GlycanTopoEqual()
    gcompeq = GlycanCompEqual()
    subsumption = GlycanSubsumption()
    topology = Topology()

    acc1 = sys.argv[1]
    acc2 = sys.argv[2]
    g1 = gtc.getGlycan(acc1)
    g2 = gtc.getGlycan(acc2)
    tg2 = topology(g2)

    verbose = True
    subshow(acc1,acc2,"<=",subsumption.leq(g1,tg2))
    subshow(acc1,acc2,"==",geq.eq(g1,tg2))

    
