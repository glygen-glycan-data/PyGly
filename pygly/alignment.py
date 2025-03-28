
from . Glycan import Glycan
from . Monosaccharide import Linkage, Substituent, Monosaccharide, Mod, Anomer, Stem
from . combinatorics import itermatchings, iterecmatchings, itergenmatchings, iterplacements, itergenmaximalmatchings, choose

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
      print("[%s, %s] %s:%s: %s"%(time.asctime(),info.function,os.path.split(info.filename)[1],info.lineno,msg))
  else:
      print("[%s, %s] %s:%s"%(time.asctime(),info.function,os.path.split(info.filename)[1],info.lineno))

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


    def parent_and_child_id_check(self, i, j, idmap=[]):

        if i.parent().is_monosaccharide() != j.parent().is_monosaccharide():
            return False

        if i.child().is_monosaccharide() != j.child().is_monosaccharide():
            return False
            
        idmapdict = dict((t[0].id(),t[1].id()) for t in idmap)

        if i.parent().is_monosaccharide():
            if idmapdict.get(i.parent().id()) != j.parent().id():
                return False
        else:
            # The monosaccharide subst attaches to.
            assert len(i.parent().parent_links()) == 1
            assert len(j.parent().parent_links()) == 1
            if idmapdict.get(i.parent().any_parent_link().parent().id()) != j.parent().any_parent_link().parent().id():
                return False

        if idmapdict.get(i.child().id()) != j.child().id():
            return False

        return True

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

    def idmap(self,a,b,leq=False):
        assert a.is_monosaccharide() == b.is_monosaccharide()
        if a.is_monosaccharide():
            if leq:
                assert self.leq(a,b)
            else:
                assert self.eq(a,b)
        else:
            assert self.substeq(a,b)
        retval = [(a.external_descriptor_id(),b.external_descriptor_id())]
        for ii,jj in itermatchings(a.substituent_links(),b.substituent_links(),
                                   lambda i,j: self.sublinkeq(i,j) and self.substeq(i.child(),j.child())):
            for i,j in zip(ii,jj): 
                if i.child().external_descriptor_id() and j.child().external_descriptor_id():
                    retval.append((i.child().external_descriptor_id(),j.child().external_descriptor_id()))
            break
        return retval

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
            la_upper = a.parent().any_parent_link()
            lb_upper = b.parent().any_parent_link()
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
            la_upper = a.parent().any_parent_link()
            lb_upper = b.parent().any_parent_link()
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
        if instonly:
            iterlinks = r.links()
        else:
            iterlinks = r.links_with_uninstantiated()
        for l in iterlinks:
            ch = l.child()
            dist.update(_mindistfromroot(ch,d+1,dist,instonly))
    return dist

class GlycanEquivalence(Comparitor):

    ### Assumes glycans have the same topology...

    def __init__(self,monocmp=None,linkcmp=None,rootmonocmp=None,substcmp=None,**kw):
        self._monocmp = monocmp
        if rootmonocmp:
            self._rootmonocmp = rootmonocmp
        else:
            self._rootmonocmp = monocmp
        self._substcmp = substcmp
        self._linkcmp = linkcmp
        self.adist = None
        self.bdist = None
        self._repeat_max_depth = 10
        self._exactmatch = True
        self._check_uninstantiated_links = True
        if "exactmatch" in kw:
            assert kw["exactmatch"] in [True, False]
            self._strict = kw["exactmatch"]
        super(GlycanEquivalence,self).__init__(**kw)

    def rootmonoeq(self,a,b):
        return self._rootmonocmp.eq(a,b)

    def monoeq(self,a,b):
        if a.is_monosaccharide() and b.is_monosaccharide():
            return self._monocmp.eq(a,b)
        elif not a.is_monosaccharide() and not b.is_monosaccharide():
            return self._substcmp.eq(a,b)
        return False

    def monoidmap(self,a,b):
        return self._monocmp.idmap(a,b)

    def linkeq(self,a,b):
        return self._linkcmp.eq(a,b)

    def repeat_max_depth(self):
        return self._repeat_max_depth

    def set_repeat_max_depth(self, d):
        assert isinstance(d, int)
        assert d > 0
        self._repeat_max_depth = d

    def subtree_eq(self,a,b,root=True,idmap=None,repeat_depth_a=0,repeat_depth_b=0):
        if self._exactmatch:
            return self.subtree_eq_strict(a,b,root=root,idmap=idmap)
        else:
            return self.subtree_eq_at_least_one(a,b,root=root,idmap=idmap,repeat_depth_a=repeat_depth_a,repeat_depth_b=repeat_depth_b)

    def subtree_eq_at_least_one(self,a,b,root=True,idmap=None,repeat_depth_a=0,repeat_depth_b=0):

        # Loose equality match
        # Return True as long as there is one chance a glycan is same as the other
        # Example:
        # -[-a-b-]- vs a-b / a-b-a-b yes
        #
        #    b
        #    |
        # -[-a-]- vs a-b yes


        if root:
            if not self.rootmonoeq(a,b):
                return False
        else:
            if not self.monoeq(a,b):
                return False

        # This assert is (and was) always trivially true, the expression evaluates to zero always 
        # - default links method did not return any repeat_bridge links
        # assert len(filter(lambda l:l.repeat_bridge_link(), a.links())) in [0, 1]
        if not a.links_has_repeat_bridge():
            alinks_no_repeat_path = a.links()
            alinks_repeat_path = []
        else:
            alinks_no_repeat_path = a.links_with_repeat_exit()
            alinks_repeat_path = a.links_with_repeat_bridge()

        # assert len(filter(lambda l:l.repeat_bridge_link(), b.links())) in [0, 1]
        if not b.links_has_repeat_bridge():
            blinks_no_repeat_path = b.links()
            blinks_repeat_path = []
        else:
            blinks_no_repeat_path = b.links_with_repeat_exit()
            blinks_repeat_path = b.links_with_repeat_bridge()

        idmaps = defaultdict(list)
        for ii,jj in itermatchings(alinks_no_repeat_path, blinks_no_repeat_path, lambda i,j: self.linkeq(i,j) and self.subtree_eq_at_least_one(i.child(),j.child(),root=False,idmap=idmaps[i,j],repeat_depth_a=repeat_depth_a,repeat_depth_b=repeat_depth_b)):
            if idmap is not None:
                idmap.append((a,b))
                for i,j in zip(ii,jj):
                    idmap.extend(idmaps[i,j])
            return True

        limit = self.repeat_max_depth()
        if repeat_depth_a > limit or repeat_depth_b > limit:
           return False

        repeat_match_a, repeat_match_b = False, False

        idmaps = defaultdict(list)
        for ii, jj in itermatchings(alinks_no_repeat_path, blinks_repeat_path,lambda i,j: self.linkeq(i,j) and self.subtree_eq_at_least_one(i.child(),j.child(),root=False,idmap=idmap,repeat_depth_a=repeat_depth_a,repeat_depth_b=repeat_depth_b+1)):
            repeat_match_b = True
            if idmap is not None:
                idmap.append((a,b))
                for i,j in zip(ii,jj):
                    idmap.extend(idmaps[i,j])
            return True

        for ii,jj in itermatchings(alinks_repeat_path, blinks_no_repeat_path,lambda i,j: self.linkeq(i,j) and self.subtree_eq_at_least_one(i.child(),j.child(),root=False,idmap=idmap,repeat_depth_a=repeat_depth_a+1,repeat_depth_b=repeat_depth_b)):
            repeat_match_a = True
            if idmap is not None:
                idmap.append((a,b))
                for i,j in zip(ii,jj):
                    idmap.extend(idmaps[i,j])
            return True

        return False

    def subtree_eq_strict(self, a, b, root=True, idmap=None):

        if root:
            if not self.rootmonoeq(a, b):
                return False
        else:
            if not self.monoeq(a, b):
                return False

        alinks = a.links()
        blinks = b.links()

        if a.links_has_repeat_bridge() or b.links_has_repeat_bridge():

            alinks_r = a.links_repeat_bridge_only()
            blinks_r = b.links_repeat_bridge_only()

            tmp = False
            for ii, jj in itermatchings(alinks_r, blinks_r, lambda i, j: self.linkeq(i, j)):
                tmp = True
                break
            if not tmp:
                return False

        idmaps = defaultdict(list)
        for ii, jj in itermatchings(alinks, blinks, lambda i, j: self.linkeq(i, j) and self.subtree_eq_strict(i.child(), j.child(), root=False, idmap=idmaps[i,j])):
            if idmap is not None:
                idmap.append((a,b))
                for i,j in zip(ii,jj):
                    idmap.extend(idmaps[i,j])
            return True

        return False

    def monosaccharide_match(self,a,b):

        if not self.monoeq(a,b):
            return False
        if self._check_uninstantiated_links:
          if a.parent_link_count() != b.parent_link_count():
            return False
        else:
          if a.parent_link_count(default=False,repeat=Linkage.NON_REPEAT,inst=Linkage.INSTANTIATED) != \
               b.parent_link_count(default=False,repeat=Linkage.NON_REPEAT,inst=Linkage.INSTANTIATED):
            return False
        if a.link_count() != b.link_count():
            return False
        if self._check_uninstantiated_links and a.link_count(default=False) != b.link_count(default=False):
            return False

        assert self.adist and self.bdist
        if self._check_uninstantiated_links and self.adist[0].get(a.id()) != self.bdist[0].get(b.id()):
            return False
        if self.adist[1].get(a.id()) != self.bdist[1].get(b.id()):
            return False

        child_links_match = False
        if self._check_uninstantiated_links:
            for ii,jj in itermatchings(a.links_with_uninstantiated(),b.links_with_uninstantiated(),
                                       lambda i,j: self.linkeq(i,j) and self.monoeq(i.child(),j.child())):
                child_links_match = True
                break
        else:
            for ii,jj in itermatchings(a.links(),b.links(),
                                       lambda i,j: self.linkeq(i,j) and self.monoeq(i.child(),j.child())):
                child_links_match = True
                break
        return child_links_match

    def eq(self, a, b, idmap=None):

        assert idmap in (None,[])
        
        self.adist = None
        self.bdist = None

        if a.has_root() and b.has_root():
            if not a.undetermined() and not b.undetermined():
                # Simple topologically determined glycan
                curidmap = []
                subtree_result = self.subtree_eq(a.root(), b.root(), idmap=curidmap)

                if not subtree_result:
                    return False

                if not self._exactmatch:
                    # Well, at least one match is found
                    if idmap is not None:
                        idmap.extend(curidmap)
                    return True

                # Exact match
                if a.repeated() and b.repeated():
                    a_rpbl = filter(lambda x: x.is_repeat_bridge(), a.all_links(include_repeat=True))
                    b_rpbl = filter(lambda x: x.is_repeat_bridge(), b.all_links(include_repeat=True))

                    for m in itermatchings(a_rpbl, b_rpbl, lambda i, j: self.linkeq(i, j) and self.parent_and_child_id_check(i, j, curidmap)):
                        if idmap is not None:
                            idmap.extend(curidmap)
                        return True
                    return False

                elif a.repeated() != b.repeated():
                    return False

                else:
                    if idmap is not None:
                        idmap.extend(curidmap)
                    return True

            if not self.subtree_eq(a.root(),b.root()):
                # non-composition, but might be undetermined toplogy. Determined part should match.
                return False

        if (a.repeated() and a.undetermined()) or (b.repeated() and b.undetermined()):
            # TODO Cannot handle cases like this...
            return False

        nodeset1 = list(a.all_nodes(subst=False,undet_subst=True))
        nodeset2 = list(b.all_nodes(subst=False,undet_subst=True))

        if len(nodeset1) != len(nodeset2):
            return False

        if a.has_root() != b.has_root():
            return False

        # if there are roots, their equivalence has been tested in subtree_align using self.rootmonoeq
        if a.has_root() and b.has_root():
            nodeset1.remove(a.root())
            nodeset2.remove(b.root())

        # we assume each node pair has a single link between it
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

        if self._check_uninstantiated_links and len(linkset1) != len(linkset2):
            return False

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
                if self._check_uninstantiated_links or not l1.instantiated():
                    continue
                l2 = linkset2.get((matching[f],matching[t]))
                if not l2 or not self.linkeq(l1,l2) or not l2.instantiated():
                    good = False
                    break
            
            if good:
                # print >>sys.stderr, "%d iterations to find an isomorphism"%(iters,)
                if idmap is not None:
                    if a.has_root() and b.has_root():
                        idmap.append((a.root(),b.root()))
                    idmap.extend(zip(ii,jj))
                self.adist = None
                self.bdist = None
                return True

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

    def monoidmap(self,a,b):
        return self._monocmp.idmap(a,b)

    def eq(self,a,b,idmap=None):

        assert idmap in (None,[])

        nodeset1 = list(a.all_nodes(subst=False,undet_subst=True))
        nodeset2 = list(b.all_nodes(subst=False,undet_subst=True))

        if len(nodeset1) != len(nodeset2):
            return False

        for ii,jj in iterecmatchings(nodeset1, nodeset2, self.nodeeq):
            # matching = dict(zip(map(lambda m: m.id(),ii),map(lambda m: m.id(),jj)))
            if idmap is not None:
                idmap.extend(zip(ii,jj))
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
        # self._topocmp = GlycanTopoEqual()

        super(GlycanPartialOrder,self).__init__(**kw)

    def compleq(self,a,b,idmap=None):
        return self._glycompcmp.leq(a,b,idmap=idmap)

    # def topoeq(self,a,b):
    #     return self._topocmp.eq(a,b)

    def rootmonoleq(self,a,b):
        return self._rootmonocmp.leq(a,b)

    def monoleq(self,a,b):
        return self._monocmp.leq(a,b)

    def monoidmap(self,a,b):
        return self._monocmp.idmap(a,b,leq=True)

    def linkleq(self,a,b):
        return self._linkcmp.leq(a,b)

    def monosaccharide_leq(self,a,b):

        if not self.monoleq(a,b):
            return False
        if a.parent_link_count(default=False,repeat=-Linkage.REPEAT_BRIDGE,inst=Linkage.INSTANTIATED) < \
             b.parent_link_count(default=False,repeat=-Linkage.REPEAT_BRIDGE,inst=Linkage.INSTANTIATED):
            return False
        if a.parent_link_count() > b.parent_link_count():
            return False
        if a.link_count() < b.link_count():
            return False
        if a.link_count(default=False,repeat=-Linkage.REPEAT_BRIDGE) > b.link_count(default=False,repeat=-Linkage.REPEAT_BRIDGE):
            return False
        assert self.adist and self.bdist
        # Using uninstantiated edges, a cannot be closer to the root
        if self.adist[0].get(a.id(),1e+20) < self.bdist[0].get(b.id(),1e+20):
            return False
        # Using instantiated edges only, b cannot be closer to the root
        if self.adist[1].get(a.id(),1e+20) > self.bdist[1].get(b.id(),1e+20):
            return False
        return True


    def subtree_leq(self, a, b, root=True, idmap=None):
        if root:
            if not self.rootmonoleq(a, b):
                return False
        else:
            if not self.monoleq(a, b):
                return False

        if idmap is not None:
            idmap.append((a,b))

        alinks = a.links()
        blinks = b.links()

        if a.links_has_repeat_bridge() or b.links_has_repeat_bridge():

            alinks_r = a.links_repeat_bridge_only()
            blinks_r = b.links_repeat_bridge_only()

            tmp = False
            for ii, jj in itermatchings(alinks_r, blinks_r, lambda i, j: self.linkeq(i, j)):
                tmp = True
            if not tmp:
                return False

        for ii, jj in itermatchings(alinks, blinks, lambda i, j: self.linkleq(i, j) and self.subtree_leq(i.child(), j.child(), root=False, idmap=idmap)):
            return True

        return False

    def leq(self,a,b,idmap=None):

        self.adist = None
        self.bdist = None

        assert idmap in (None,[])

        lineno("Start of Subsumption leq")

        # a's composition must be subsumed by b's composition
        curidmap = []
        if not self.compleq(a,b,idmap=curidmap):
            return False

        lineno("Composition is subsumed")

        # if b is a composition, it doesn't matter what a's topology is
        if not b.has_root():
            if idmap is not None:
                idmap.extend(curidmap)
            return True

        lineno("b is not a composition")

        # b is not a composition, so a must also not be a composition
        if not a.has_root():
            return False

        # Both a and b are now rooted topologies, though not nec. topologically determined
        
        lineno()

        if a.repeated() != b.repeated():
            return False

        if not a.undetermined() and not b.undetermined():
            # Simple topologically determined glycan
            curidmap = []
            subtree_result = self.subtree_leq(a.root(),b.root(), idmap=curidmap)

            if not subtree_result:
                return False

            if a.repeated() and b.repeated():
                a_rpbl = filter(lambda x: x.is_repeat_bridge(), a.all_links(include_repeat=True))
                b_rpbl = filter(lambda x: x.is_repeat_bridge(), b.all_links(include_repeat=True))

                for m in itermatchings(a_rpbl, b_rpbl, lambda i, j: self.linkleq(i,j) and self.parent_and_child_id_check(i, j, curidmap)):
                    if idmap is not None:
                        idmap.extend(curidmap)
                    return True
                return False

            elif a.repeated() != b.repeated():
                return False

            else:
                if idmap is not None:
                    idmap.extend(curidmap)
                return True

        if (a.repeated() and a.undetermined()) or (b.repeated() and b.undetermined()):
            # TODO Cannot handle cases like this...
            return False

        lenns1uds = len(list(filter(lambda n: not n.is_monosaccharide(),a.all_nodes(subst=False,undet_subst=True))))
        lenns2uds = len(list(filter(lambda n: not n.is_monosaccharide(),b.all_nodes(subst=False,undet_subst=True))))
        if lenns2uds > lenns1uds:
            a = a.clone()
        elif lenns2uds < lenns1uds:
            return False

        self.adist = _mindistsfromroot(a.root())
        self.bdist = _mindistsfromroot(b.root())

        lineno()

        nodeset1 = list(a.all_nodes(subst=False))
        nodeset2 = list(b.all_nodes(subst=False))

        nodeset1all = list(a.all_nodes(subst=False,undet_subst=True))
        nodeset2all = list(b.all_nodes(subst=False,undet_subst=True))

        nodeset1uds = list(filter(lambda n: not n.is_monosaccharide(),nodeset1all))
        nodeset2uds = list(filter(lambda n: not n.is_monosaccharide(),nodeset2all))

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
                key = (l.parent().any_parent_link().parent().id(), l.child().id())
            assert (key not in linkset1)
            linkset1[key] = l

        linkset2 = dict()
        for l in b.all_links(uninstantiated=True):
            if l.parent().is_monosaccharide():
                key = (l.parent().id(), l.child().id())
            else:
                key = (l.parent().any_parent_link().parent().id(), l.child().id())
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

        # print(allsubst)
        # print(subst1link)
        # print(subst1float)
        # print(subst2link)
        # print(subst2float)

        subst2move = []
        bad = False
        for s in allsubst:
            if (subst1link[s]+subst1float[s]) != (subst2link[s]+subst2float[s]):
                return False
            elif subst2float[s] > subst1float[s] and \
                     subst2link[s] == 0:
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
            raise RuntimeError("Cannot resolve partial order, floating substituents too complicated")

        # print subst2move
            
        # strip specific floating substituents from nodes in nodeset1
        newnodeset1 = []
        for m1 in nodeset1:
            for sl in list(m1.substituent_links()):
                if sl.child().name() in subst2move:
                    m1.remove_substituent_link(sl)
                    sl.child().del_parent_link(sl)
            newnodeset1.append(m1)

            # print newnodeset1

        lineno()

        iters = 0
        for ii,jj in itergenmatchings(newnodeset1, nodeset2, self.monosaccharide_leq):

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
                    if (f,t) not in matchedlinkset2links and linkset2[(f,t)].instantiated():
                        good = False
                        break
                
            if good:
                # print >>sys.stderr, "%d iterations to find an isomorphism"%(iters,)
                if idmap is not None:
                    idmap.append((a.root(),b.root()))
                    idmap.extend(zip(ii,jj))
                return True
                
        lineno()

        # print >>sys.stderr, "%d iterations to prove no isomorphism exists"%(iters,)
        return False



class ConnectedNodesCache:

    def __init__(self):
        self.clear()

    def put(self, g):
        self.clear()
        for m in g.all_nodes():
            self.data[m] = defaultdict(list)
            self.data[m][1].append({m})
            for c in set(l.child() for l in m.links_with_uninstantiated()):
                self.data[m][2].append({m,c})

    def update_cache(self, m, size):

        if size in self.data[m]:
            return

        i = max(self.data[m].keys())

        while i < size:
            i += 1
            res = dict()
            for currentSet in self.data[m][i-1]:
                for res0 in self.connectedNodesPlusOneSimple(currentSet):
                    res0sig = frozenset(res0)
                    if res0sig not in res:
                        res[res0sig] = res0
            self.data[m][i] = list(res.values())

    def connectedNodesPlusOneSimple(self, currentSet):
        res = set(x.child() for n in currentSet for x in n.links_with_uninstantiated())
        children = res - currentSet

        res = []
        for n in children:
            res.append(currentSet|{n})
        return res

    def get(self, m, size):
        size = int(size)
        self.update_cache(m, size)
        return self.data[m][size]

    def clear(self):
        self.data = {}

class SubstructureSearch(GlycanPartialOrder):

    def __init__(self, **kw):
        self.connected_nodes_pre_computed = kw.get("connected_nodes_pre_computed", True)
        if self.connected_nodes_pre_computed:
            self.nodes_cache = kw.get("connected_nodes_cache", ConnectedNodesCache())
        self._repeat_max_depth = 10
        self._glyeq = GlycanEqual()
        super(SubstructureSearch, self).__init__(**kw)

    def repeat_max_depth(self):
        return self._repeat_max_depth

    def set_repeat_max_depth(self, d):
        assert isinstance(d, int)
        assert d > 0
        self._repeat_max_depth = d

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
                for ll in tgp.links_with_uninstantiated():
                    if ll.child() not in tg_nodes:
                        continue

                    if motif_nodes.index(c) == tg_nodes.index(ll.child()) and self._linkcmp.leq(l, ll):

                        if not self.check_links_childlink_count(c, ll.child()):
                            # No real usage in class SubstructureSearch
                            # Special case for Non-reducing end matcher
                            continue

                        found_matched_link = True
                        break

                if not found_matched_link:
                    return False

        return True

    def check_links_childlink_count(self, m, tg):
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
            for cl in c.links_with_uninstantiated():
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
                                      set([x.child() for x in r.links_with_uninstantiated()]),
                                      size):
            res.append(cn)

        return res



    def allConnectedNodesByRoot(self, r, size):
        if self.connected_nodes_pre_computed:
            return self.nodes_cache.get(r, size)
        else:
            return self.allConnectedNodesByRootRecursive(r, size)

    def subtree_leq(self, m, tg, root=True, repeat_depth_m=0, repeat_depth_tg=0, idmaps=None):

        if root:
            if not self.rootmonoleq(m, tg):
                return False
        else:
            if not self.monoleq(m, tg):
                return False

        # print m.id(),tg.id()

        ######
        # *_no_repeat_edges includes NON_REPEAT and REPEAT_EXIT
        # *_repeat_edges includes NON_REPEAT AND REPEAT_BRIDGE
        # *_nonrepeat_edges includes NON_REPEAT only
        # *_repeat_bridge_only_edges includes REPEAT_BRIDGE only
        #
        # NOTE: Assumes at most one repeat bridge edge for each
        # monosaccharide and never aligns repeat bridge to repeat
        # bridge (should it?)  Having a repeat bridge for motif and
        # glycan monosaccharide *should* be rare!!!!
        # 
        # Asymmetry in motif and glycan edge combinatorics is because
        # we dont have to use all fo the glycan monosaccharide's
        # edges, but we must use all the motif monosaccharide's edges
        ######

        assert len(list(m.links_repeat_bridge_only())) <= 1
        if not m.links_has_repeat_bridge():
            mlinks_no_repeat_edges = list(m.links()) # same as links_with_repeat_exit()
            mlinks_repeat_edges = []
        else:
            mlinks_no_repeat_edges = list(m.links_with_repeat_exit())
            mlinks_repeat_edges = list(m.links_with_repeat_bridge())

        assert len(list(tg.links_repeat_bridge_only())) <= 1
        if not tg.links_has_repeat_bridge():
            tglinks_no_repeat_edges = list(tg.links()) # same as links_with_repeat_exit()
            tglinks_repeat_bridge_only_edges = []
        else:
            tglinks_nonrepeat_edges = list(tg.links_without_repeat())
            tglinks_no_repeat_edges = list(tg.links_with_repeat_exit())
            tglinks_repeat_bridge_only_edges = list(tg.links_repeat_bridge_only())

        for tglinks in choose(tglinks_no_repeat_edges, len(mlinks_no_repeat_edges)):
            recidmaps = defaultdict(list)
            for ii,jj in itermatchings(mlinks_no_repeat_edges, tglinks, lambda i,j: self.linkleq(i,j) and self.subtree_leq(i.child(),j.child(),root=False,repeat_depth_m=repeat_depth_m,repeat_depth_tg=repeat_depth_tg,idmaps=recidmaps[i,j])):
                if idmaps is not None:
                    thisidmaps = [[(m,tg)]]
                    for i,j in zip(ii,jj):
                        idmaps1 = []
                        for idmap0 in thisidmaps:
                            for idmap1 in recidmaps[i,j]:
                                idmaps1.append(idmap0 + idmap1)
                        thisidmaps = idmaps1
                    idmaps.extend(thisidmaps)
                else:
                    return True

        # print "HERE!0"
        # print m.id(),tg.id()
        # print [",".join(map(lambda t: str(t[0].id())+":"+str(t[1].id()),idmap)) for idmap in idmaps]

        limit = self.repeat_max_depth()
        if repeat_depth_m >= limit or repeat_depth_tg >= limit:
            if len(idmaps) > 0:
                return True
            return False

        # print "HERE!1"
        # print m.id(),tg.id()
        # print [",".join(map(lambda t: str(t[0].id())+":"+str(t[1].id()),idmap)) for idmap in idmaps]

        if len(mlinks_no_repeat_edges) > 0 and len(tglinks_repeat_bridge_only_edges) > 0:
            for tglinks in choose(tglinks_nonrepeat_edges, len(mlinks_no_repeat_edges)-1):
                tglinks += tglinks_repeat_bridge_only_edges
                recidmaps = defaultdict(list)
                for ii, jj in itermatchings(mlinks_no_repeat_edges, tglinks, lambda i,j: self.linkleq(i,j) and self.subtree_leq(i.child(),j.child(),root=False,repeat_depth_m=repeat_depth_m,repeat_depth_tg=repeat_depth_tg+1,idmaps=recidmaps[i,j])):
                    if idmaps is not None:
                        thisidmaps = [[(m,tg)]]
                        for i,j in zip(ii,jj):
                            idmaps1 = []
                            for idmap0 in thisidmaps:
                                for idmap1 in recidmaps[i,j]:
                                    idmaps1.append(idmap0 + idmap1)
                            thisidmaps = idmaps1
                        idmaps.extend(thisidmaps)
                    else:
                        return True

        # print "HERE!2"
        # print m.id(),tg.id()
        # print [",".join(map(lambda t: str(t[0].id())+":"+str(t[1].id()),idmap)) for idmap in idmaps]

        if len(mlinks_repeat_edges) > 0 and len(tglinks_no_repeat_edges) > 0:
          for tglinks in choose(tglinks_no_repeat_edges, len(mlinks_repeat_edges)):
            recidmaps = defaultdict(list)
            for ii,jj in itermatchings(mlinks_repeat_edges, tglinks,lambda i,j: self.linkleq(i,j) and self.subtree_leq(i.child(),j.child(),root=False,repeat_depth_m=repeat_depth_m+1,repeat_depth_tg=repeat_depth_tg,idmaps=recidmaps[i,j])):
                if idmaps is not None:
                    thisidmaps = [[(m,tg)]]
                    for i,j in zip(ii,jj):
                        idmaps1 = []
                        for idmap0 in thisidmaps:
                            for idmap1 in recidmaps[i,j]:
                                idmaps1.append(idmap0 + idmap1)
                        thisidmaps = idmaps1
                    idmaps.extend(thisidmaps)
                else:
                    return True


        # print m.id(),tg.id()
        # print [",".join(map(lambda t: str(t[0].id())+":"+str(t[1].id()),idmap)) for idmap in idmaps]

        if len(idmaps) > 0:
            return True
        return False


    def leq(self, m, tg, rootOnly=False, anywhereExceptRoot=False, underterminedLinkage=True, idmaps=None):

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

        if idmaps is not None:
            seenidmaps = set()

        # Use subtree algorithm to save runtime
        if not m.undetermined():
            for n in potential_TG_root:
                curidmaps = []
                if self.subtree_leq(m.root(), n, idmaps=curidmaps):
                    if idmaps is None:
                        return True
                    for idmap in curidmaps:
                        idmapkey = tuple(sorted(set(map(lambda t: (t[0].id(),t[1].id()),idmap))))
                        if idmapkey not in seenidmaps:
                            idmaps.append(idmap)
                            seenidmaps.add(idmapkey)
            
            # return False when no match based on subtree algorithm and the glycan is not undetermined.
            if not tg.undetermined():
                if idmaps is None:
                    return False
                return len(idmaps)>0

        if not underterminedLinkage:
            if idmaps is None:
                return False
            return len(idmaps)>0

        if tg.undetermined() and self.connected_nodes_pre_computed:
            # pre compute the connected node set
            self.nodes_cache.put(tg)

        undetroots = [ sorted(ec,key=lambda m: m.id()) for ec in tg._undetermined ]
        undetidmap = dict()
        for tgm in tg.all_nodes():
            undetidmap[tgm] = tgm.id()
        for ec in undetroots:
            for udr in ec[1:]:
                if udr.is_monosaccharide():
                    idmap = []
                    if self._glyeq.eq(Glycan(ec[0]),Glycan(udr),idmap=idmap):
                        for m0,m1 in idmap:
                            undetidmap[m1] = m0.id()
                else:
                    undetidmap[udr] = ec[0].id()

        #i = 0
        for tg_root in potential_TG_root:

            skip = 0
            check = 0
            generate = 0
            tgsetsigs = dict()
            for tg_nodes in self.allConnectedNodesByRoot(tg_root, len(list(m.all_nodes()))):
                tgsetsig = tuple(sorted(map(undetidmap.get,tg_nodes)))
                if tgsetsigs.get(tgsetsig) is not None:
                    if tgsetsigs.get(tgsetsig) != False:
                        for idmapi in tgsetsigs[tgsetsig]:
                            tgn0=[ t[1] for t in idmapi ]
                            # print("a",sorted([ t.id() for t in tgn0 ]))
                            # print("b",sorted([ t.id() for t in tg_nodes ]))
                            for ii,jj in itergenmatchings(tgn0,tg_nodes,
                                              lambda a,b: undetidmap[a] == undetidmap[b]):
                                if ii == jj:
                                    continue
                                tgnmap = dict(zip(ii,jj))
                                idmap1 = [ (t[0],tgnmap[t[1]]) for t in idmapi ]
                                idmapkey = tuple(sorted(set(map(lambda t: (t[0].id(),t[1].id()),idmap1))))
                                if idmapkey not in seenidmaps:
                                    idmaps.append(idmap1)
                                    seenidmaps.add(idmapkey)
                        generate += 1
                    else:
                        skip += 1
                    continue

                check += 1
                tgsetsigs[tgsetsig] = False
                # print(check,skip,generate)

                for monoset_m, monoset_tg in itergenmatchings(list(m.all_nodes()), tg_nodes, self.monoleq):

                    if self.check_links(m, monoset_m, monoset_tg):
                        if idmaps is None:
                            return True
                        idmap = list(zip(monoset_m, monoset_tg))
                        idmapkey = tuple(sorted(set(map(lambda t: (t[0].id(),t[1].id()),idmap))))
                        if idmapkey not in seenidmaps:
                            idmaps.append(idmap)
                            seenidmaps.add(idmapkey)
                            if tgsetsigs[tgsetsig] == False:
                                tgsetsigs[tgsetsig] = [ idmap ]
                            else:
                                tgsetsigs[tgsetsig].append(idmap)

        if idmaps is not None:
            return len(idmaps)>0

        return False

    def whole_glycan_match_check(self, m, tg):

        if m.repeated() or tg.repeated():
            return False

        if len(list(m.all_nodes())) != len(list(tg.all_nodes())):
            return False

        if len(list(m.all_nodes(subst=True))) != len(list(tg.all_nodes(subst=True))):
            return False

        if m.underivitized_molecular_weight() != tg.underivitized_molecular_weight():
            return False

        return True

    def whole_glycan_match(self, m, tg):

        if not self.whole_glycan_match_check(m, tg):
            return True

        return self.leq(m, tg, rootOnly=True)

    def idmaps_toids(self,idmaps):
        newidmaps = []
        for idmap in idmaps:
            idmapids = [ ti for t in idmap for ti in self.monoidmap(*t) ]
            newidmaps.append(idmapids)
        return newidmaps

    def matched_ids(self,idmaps,glycan):
        allstructids = set()
        allstructlinkids = set()
        for idmap in idmaps: 
            structids = set(t[1] for t in idmap)
            allstructids.update(structids)
            monoids = set(t1 for t1 in structids if '.' not in t1)
            linkids = set()
            if glycan:
                for l in glycan.all_links(uninstantiated=True): 
                    if l.parent().external_descriptor_id() in monoids and \
                       l.child().external_descriptor_id() in monoids:
                        if l.instantiated():
                            linkids.add((l.parent().external_descriptor_id(),l.child().external_descriptor_id())) 
                        else:
                            linkids.add(('',l.child().external_descriptor_id())) 
                # if len(monoids) != (len(linkids) + 1): 
                #     print("Warning: Bad linkids length motifacc: %s structacc: %s"%(motifacc,structacc),file=sys.stderr) 
                allstructlinkids.update(linkids) 
        return sorted(str(i) for i in allstructids), sorted("%s-%s"%p for p in allstructlinkids)

class SubstructureSearchNonReducingEnd(SubstructureSearch):

    def check_links_childlink_count(self, m, tg):

        # Reject condition
        if not m.has_links() and tg.has_links():
            return False

        return True

    def subtree_leq(self, m, tg, root=True, **karg):

        if not self.check_links_childlink_count(m, tg):
            return False

        return super(SubstructureSearchNonReducingEnd, self).subtree_leq(m, tg, root=root, **karg)




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

    def leq(self,a,b,idmap=None):
      
        assert idmap in (None,[])

        lineno()

        nodeset1 = list(a.all_nodes(subst=False))
        nodeset2 = list(b.all_nodes(subst=False))

        nodeset1all = list(a.all_nodes(subst=False,undet_subst=True))
        nodeset2all = list(b.all_nodes(subst=False,undet_subst=True))

        nodeset1uds = list(filter(lambda n: not n.is_monosaccharide(),nodeset1all))
        nodeset2uds = list(filter(lambda n: not n.is_monosaccharide(),nodeset2all))

        lineno()

        if len(nodeset1) != len(nodeset2):
            return False

        lineno()

        # print(len(nodeset1),len(nodeset1uds),len(nodeset2),len(nodeset2uds))

        # for m in nodeset1all:
        #     print(m)

        # for m in nodeset2all:
        #     print(m)

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
                if idmap is not None:
                    idmap.extend(zip(ii,jj))
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

            # print(allsubst)
            # print(subst1link)
            # print(subst1float)
            # print(subst2link)
            # print(subst2float)

            subst2move = []
            bad = False
            for s in allsubst:
                if (subst1link[s]+subst1float[s]) != (subst2link[s]+subst2float[s]):
                    return False
                elif subst2float[s] > subst1float[s] and \
                     subst2link[s] == 0:
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
                for sl in list(m1.substituent_links()):
                    if sl.child().name() in subst2move:
                        m1.remove_substituent_link(sl)
                        sl.child().del_parent_link(sl)
                newnodeset1.append(m1)

            # print newnodeset1
            # print nodeset2

            # Now the monosaccharides should match....
            for ii,jj in itergenmatchings(newnodeset1,nodeset2,self.nodeleq):
                if idmap is not None:
                    idmap.extend(zip(ii,jj))
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
    # Strict

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


class MonosaccharideMotifComparisonGlyTouCanOriginal(MonosaccharideComparitor):

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


class MonosaccharideMotifComparisonSubstTolerance(MonosaccharideComparitor):

    def leq(self, m, g):
        if m._anomer and g._anomer and m._anomer != g._anomer:
            return False
        if m._config and m._config != g._config:
            return False
        if m._stem and m._stem != g._stem:
            return False
        if m._superclass != g._superclass:
            return False
        if m._ring_start and g._ring_start and m._ring_start != g._ring_start:
            return False
        if m._ring_end and g._ring_end and m._ring_end != g._ring_end:
            return False

        gmod = copy.deepcopy(g._mods)
        for mod in m._mods:
            if mod in gmod:
                gmod.remove(mod)
            else:
                return False

        for mod in gmod:
            # aldi tolerance
            if mod[1] != Mod.aldi:
                return False

        any = False

        msl = list(m.substituent_links())
        gsl = list(g.substituent_links())
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

    def idmap(self,a,b,**kw):
        assert a.is_monosaccharide() == b.is_monosaccharide()
        if a.is_monosaccharide():
            assert self.leq(a,b)
        else:
            assert self.substeq(a,b)

        retval = [(a.external_descriptor_id(),b.external_descriptor_id())]
        any = False

        asl = list(a.substituent_links())
        bsl = list(b.substituent_links())
        if len(asl) > len(bsl):
            return False

        bsl_mandatory, bsl_optional = [], []
        for sl in bsl:
            if sl.child().name() in [Substituent.sulfate, Substituent.phosphate]:
                bsl_optional.append(sl)
            else:
                bsl_mandatory.append(sl)

        for x in choose(bsl_optional, len(asl) - len(bsl_mandatory)):

            for ii, jj in itermatchings(asl, x + bsl_mandatory,
                                        lambda i, j: self.sublinkeq(i, j) and self.substeq(i.child(), j.child())):
                for i,j in zip(ii,jj): 
                    if i.child().external_descriptor_id() and j.child().external_descriptor_id():
                        retval.append((i.child().external_descriptor_id(),j.child().external_descriptor_id()))
                any = True
                break
            if any:
                break
        assert(any)
        return retval

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
        if a._link_type != b._link_type:
            return False
        return True

class LinkageTopoEqualSimple(LinkageComparitorBase):
    def eq(self,a,b):
        if a._parent_type != b._parent_type:
            return False
        if a._child_type != b._child_type:
            return False
        # do not check parent_pos or child_pos
        if a._link_type != b._link_type:
            return False
        return True

class LinkageImageEqualSimple(LinkageComparitorBase):
    def eq(self,a,b):
        if a._link_type != b._link_type:
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
        if not a.instantiated() and b.instantiated():
            return False
        return True

class LinkageTopoSubsumedSimple(LinkageSubsumedSimple):

    def leq(self,a,b):
        if not self._leq_(a._parent_type,b._parent_type):
            return False
        if not self._leq_(a._child_type,b._child_type):
            return False
        if not a.instantiated() and b.instantiated():
            return False
        return True

class LinkageComparitorMotifGlyTouCanOriginalSimple(LinkageComparitorBase):

    def leq(self, m, g):
        if m._parent_type != g._parent_type:
            return False
        if m._child_type != g._child_type:
            return False

        if m._child_pos and g._child_pos and m._child_pos.isdisjoint(g._child_pos):
            return False

        if m._parent_pos and g._parent_pos and m._parent_pos.isdisjoint(g._parent_pos):
            return False

        return True

class LinkageComparitorMotifStrictSimple(LinkageComparitorBase):

    def leq(self, m, g):
        if m._parent_type != g._parent_type:
            return False
        if m._child_type != g._child_type:
            return False

        if m._child_pos and not g._child_pos:
            return False
        elif m._child_pos and g._child_pos and not (m._child_pos >= g._child_pos):
            return False

        if m._parent_pos and not g._parent_pos:
            return False
        elif m._parent_pos and g._parent_pos and not (m._parent_pos >= g._parent_pos):
            return False

        return True

class SubstLinkageComparitorMotifSimple(LinkageEqualSimple):

    def leq(self,a,b):
        return self.eq(a,b)

class SubstLinkageComparitorMotifStrictSimple(LinkageComparitorMotifStrictSimple):

    def eq(self,a,b):
        return self.leq(a,b)

class SubstLinkageComparitorMotifLooseSimple(LinkageComparitorMotifGlyTouCanOriginalSimple):

    def eq(self,a,b):
        return self.leq(a,b)

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

class LinkageComparitorMotifGlyTouCanOriginal(LinkageComparitor):

    def __init__(self, **kw):
        kw['substcmp'] = SubstituentEqual(**kw)
        kw['linkcmp'] = LinkageComparitorMotifGlyTouCanOriginalSimple(**kw)
        kw['sublinkcmp'] = SubstLinkageComparitorMotifSimple(**kw)
        super(LinkageComparitorMotifGlyTouCanOriginal, self).__init__(**kw)

class LinkageComparitorMotifStrict(LinkageComparitor):

    def __init__(self, **kw):
        kw['substcmp'] = SubstituentEqual(**kw)
        kw['linkcmp'] = LinkageComparitorMotifStrictSimple(**kw)
        kw['sublinkcmp'] = SubstLinkageComparitorMotifStrictSimple(**kw)
        super(LinkageComparitorMotifStrict, self).__init__(**kw)


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
        self._check_uninstantiated_links = False

class GlycanCompImageEqual(CompositionEquivalence):
    def __init__(self,**kw):
        kw['substcmp']=SubstituentEqual(**kw)
        kw['sublinkcmp']=LinkageImageEqual(**kw)
        # monotest needs a subst test and a sublink test
        kw['monocmp']=MonosaccharideImageEqual(**kw)
        super(GlycanCompImageEqual,self).__init__(**kw)

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

class GlyTouCanMotifOriginal(SubstructureSearch):

    def __init__(self, **kw):
        kw["monocmp"] = MonosaccharideMotifComparisonGlyTouCanOriginal(
            substcmp=SubstituentEqual(),
            sublinkcmp=LinkageEqual()
        )
        kw["linkcmp"] = LinkageComparitorMotifGlyTouCanOriginal()
        super(GlyTouCanMotifOriginal, self).__init__(**kw)


class GlyTouCanMotifNonReducingEndOriginal(SubstructureSearchNonReducingEnd):

    def __init__(self, **kw):
        kw["monocmp"] = MonosaccharideMotifComparisonGlyTouCanOriginal(
            substcmp=SubstituentEqual(),
            sublinkcmp=LinkageEqual()
        )
        kw["linkcmp"] = LinkageComparitorMotifGlyTouCanOriginal()
        super(GlyTouCanMotifNonReducingEndOriginal, self).__init__(**kw)

class MotifStrict(SubstructureSearch):

    def __init__(self, **kw):
        kw["monocmp"] = MonosaccharideMotifComparison(
            substcmp=SubstituentEqual(),
            sublinkcmp=SubstLinkageComparitorMotifStrictSimple()
        )
        kw["linkcmp"] = LinkageComparitorMotifStrict()
        super(MotifStrict, self).__init__(**kw)

class NonReducingEndMotifStrict(SubstructureSearchNonReducingEnd):

    def __init__(self, **kw):
        kw["monocmp"] = MonosaccharideMotifComparison(
            substcmp=SubstituentEqual(),
            sublinkcmp=SubstLinkageComparitorMotifStrictSimple()
        )
        kw["linkcmp"] = LinkageComparitorMotifStrict()
        super(NonReducingEndMotifStrict, self).__init__(**kw)

class WholeGlycanEqualMotifStrict(GlycanPartialOrder):

    def __init__(self,**kw):
        kw['substcmp']=SubstituentEqual(**kw)
        kw['linkcmp']=LinkageComparitorMotifStrict(**kw)
        kw['sublinkcmp']=SubstLinkageComparitorMotifStrictSimple(**kw)
        # monotest needs a subst test and a sublink test
        kw['monocmp']=MonosaccharideMotifComparison(**kw)
        super(WholeGlycanEqualMotifStrict,self).__init__(**kw)

class MotifInclusive(SubstructureSearch):

    def __init__(self, **kw):
        kw["monocmp"] = MonosaccharideMotifComparisonSubstTolerance(
            substcmp=SubstituentEqual(),
            sublinkcmp=SubstLinkageComparitorMotifLooseSimple()
        )
        kw["linkcmp"] = LinkageComparitorMotifGlyTouCanOriginal()
        super(MotifInclusive, self).__init__(**kw)

class NonReducingEndMotifInclusive(SubstructureSearchNonReducingEnd):
    def __init__(self, **kw):
        kw["monocmp"] = MonosaccharideMotifComparisonSubstTolerance(
            substcmp=SubstituentEqual(),
            sublinkcmp=SubstLinkageComparitorMotifLooseSimple()
        )
        kw["linkcmp"] = LinkageComparitorMotifGlyTouCanOriginal()
        super(NonReducingEndMotifInclusive, self).__init__(**kw)

class WholeGlycanEqualMotifInclusive(GlycanPartialOrder):

    def __init__(self,**kw):
        kw['substcmp']=SubstituentEqual(**kw)
        kw['linkcmp']=LinkageComparitorMotifGlyTouCanOriginal(**kw)
        kw['sublinkcmp']=SubstLinkageComparitorMotifLooseSimple(**kw)
        # monotest needs a subst test and a sublink test
        kw['monocmp']=MonosaccharideMotifComparisonSubstTolerance(**kw)
        super(WholeGlycanEqualMotifInclusive,self).__init__(**kw)

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
        print("g1(%s), g2(%s): %s = %s (%s expected): Failed!"%(acc1,acc2,test,result,expected))
        # sys.exit(1)
    else:
        print("g1(%s), g2(%s): %s = %s (%s expected)"%(acc1,acc2,test,result,expected))

def subshow(acc1,acc2,test,result):

    if result:
        print("%s %s %s"%(acc1,test,acc2))
    else:
        print("%s %s/%s %s"%(acc1,test[0],test[-1],acc2))

def compare(gd1,gd2):

    acc1,g1,tacc1,cacc1,bcacc1,l1 = list(map(gd1.get,("acc","glycan","topo","comp","bcomp","level")))
    acc2,g2,tacc2,cacc2,bcacc2,l2 = list(map(gd2.get,("acc","glycan","topo","comp","bcomp","level")))

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

    acc1,g1,tacc1,cacc1,bcacc1,l1 = list(map(gd1.get,("acc","glycan","topo","comp","bcomp","level")))
    acc2,g2,tacc2,cacc2,bcacc2,l2 = list(map(gd2.get,("acc","glycan","topo","comp","bcomp","level")))

    subshow(acc1,acc2,"g1 <= g2",subsumption.leq(g1,g2))


def main():       
    from collections import defaultdict
    from . manipulation import Topology, Composition
    global verbose

    # from GlyTouCan import GlyTouCan
    from . GlycanResource import GlyTouCan
    gtc = GlyTouCan(usecache=False,prefetch=False)

    from . GlycanFormatter import GlycoCTFormat, WURCS20Format, GlycanParseError
    from . GlycanFormatter import IUPACLinearFormat, IUPACParserExtended1
    from . CompositionFormatter import CompositionFormat
    glycoct_format = GlycoCTFormat()
    wurcs_format = WURCS20Format()
    iupac_format = IUPACLinearFormat()

    geq = GlycanEqual()
    gtopoeq = GlycanTopoEqual()
    gcompeq = GlycanCompEqual()
    subsumption = GlycanSubsumption()
    topology = Topology()

    acc1 = sys.argv[1]
    g1 = None
    if os.path.exists(acc1):
        seq = open(acc1).read()
        for fmt in (glycoct_format,wurcs_format,iupac_format):
            try:
                g1 = fmt.toGlycan(seq)
            except GlycanParseError:
                pass
    else:
        g1 = gtc.getGlycan(acc1)
    assert g1 != None

    acc2 = sys.argv[2]
    g2 = None
    if os.path.exists(acc2):
        seq = open(acc2).read()
        for fmt in (glycoct_format,wurcs_format,iupac_format):
            try:
                g2 = fmt.toGlycan(seq)
            except GlycanParseError:
                pass
    else:
        g2 = gtc.getGlycan(acc2)
    assert g2 != None

    # print(g1)
    # print(g2)
    # tg2 = topology(g2)

    verbose = True
    subshow(acc1,acc2,"<=",subsumption.leq(g1,g2))
    # subshow(acc1,acc2,"==",geq.eq(g1,tg2))
    
if __name__ == "__main__":
    main()
