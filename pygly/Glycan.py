#Glycan
#Kevin B Chandler
import operator
import sys
import time
import copy
from collections import defaultdict

try:
    from itertools import permutations, product
except ImportError:
    from combinatorics import permutations, product

from combinatorics import itermatchings, iterecmatchings

from Monosaccharide import Monosaccharide, Linkage
from MonoFormatter import IUPACSym, LinCodeSym

iupacSym = IUPACSym()
lcSym = LinCodeSym()

from CompositionTable import Composition,ResidueCompositionTable,PermethylCompositionTable
from ElementMass import MonoisotopicElementMass
from MonoFactory import MonoFactory
from MonoFormatter import MassSym

ctable = ResidueCompositionTable()
pctable = PermethylCompositionTable()
elmt = MonoisotopicElementMass()
mfactory = MonoFactory()
msym = MassSym()

class Glycan:

    iupacSym = IUPACSym()
    lcSym = LinCodeSym()

    def __init__(self,root=None):
        self.set_root(root)
        self._undetermined = None
        self._bions = None
        self._yions = None

    def root(self):
        return self._root

    def set_root(self, r):
        self._root = r

    def set_ids(self):
	for i,m in enumerate(self.all_nodes(subst=True)):
	    m.set_id(i+1)

    def unset_ids(self):
	for m in self.all_nodes(subst=True):
	    m.unset_id()

    def set_undetermined(self, und):
        if und == None or len(und) == 0:
            self._undetermined = None
            return
        u = list(und)
        ueq = defaultdict(set)
        placed = set()
        for i in range(len(u)):
            if i in placed:
                continue
            placed.add(i)
            ueq[i].add(u[i])
            for j in range(i+1,len(u)):
                if j in placed:
                    continue
		if not self.undetroot_equals(u[i],u[j],mapids=False):
		    continue
                ueq[i].add(u[j])
                placed.add(j)
        self._undetermined = list(ueq.values())

    def undetermined(self):
	return self._undetermined != None

    def undetermined_roots(self):
        if self._undetermined != None:
            for ec in self._undetermined:
                for r in ec:
                    yield r
    
    def undetermined_root_reprs(self):
        if self._undetermined != None:
            for ec in self._undetermined:
                for r in ec:
                    yield (r,len(ec))
                    break

    def unconnected_roots(self):
        for r in self.undetermined_roots():
            if not r.connected():
                yield r

    def isolated_nodes(self):
        for r in self.unconnected_roots():
            if len(r.parent_links()) == 0:
                yield r

    def isolated_node_reprs(self):
        if self._undetermined != None:
            for ec in self._undetermined:
                count = 0
                repr = None
                for r in ec:
                    if not r.connected() and len(r.parent_links()) == 0:
                        count += 1
                        if not repr:
                            repr = r
                yield (repr,count)

    def has_root(self):
        return (self._root != None)

    def fully_determined(self):
        if self.undetermined():
            return False
        for m in self.all_nodes(subst=True):
            if not m.fully_determined():
                return False
        for l in self.all_links(subst=True):
            if not m.fully_determined():
                return False
        return True

##     def add_instantiation(self, inst):
## 	if self._instantiations == None:
## 	    self._instantiations = []
## 	self._instantiations.append(inst)

##     maxlinks = {'Fuc': 0, 'NeuAc': 0, 'NeuGc': 0, 'Xyl': 0}
##     def auto_instantiations(self):
## 	undetsets = defaultdict(set)
## 	todo = [self.root()]
##         while len(todo) > 0:
##             m = todo.pop(0)
##             for l in m.substituent_links(False):
##                 if l.undetermined():
## 		    undetsets[l.child()].add(l)
##             for l in m.links(False):
##                 if l.undetermined():
## 		    undetsets[l.child()].add(l)
##                 todo.insert(0,l.child())
## 	# Pick one from each child-set
## 	for inst in product(*(undetsets.values())):
## 	    # Potentially, eliminate infeasible combinations of
## 	    # instantiated edges, too many on a parent, bond already
## 	    # used, etc.
## 	    counts = defaultdict(int)
## 	    counts1 = defaultdict(int)
## 	    for l in inst:
## 		if l.parent_pos():
## 		    counts[(l.parent(),l.parent_pos())] += 1
## 		counts1[l.parent()] += 1
## 	    for p in counts1:
## 		for l in p.links():
## 		    if l.undetermined():
## 			continue
## 		    if l.parent_pos():
## 		        counts[(l.parent(),l.parent_pos())] += 1
## 		    counts1[l.parent()] += 1
## 	    coremannose = set()
## 	    for m in self.root().children():
## 		for m1 in m.children():
## 		    try:
## 	                if iupacSym.toStr(m1) == 'Man':
## 			    coremannose.add(m1)
## 		    except KeyError:
## 			pass
## 	    # print counts
## 	    bad = False
## 	    for m,c in counts1.items():
## 		try:
## 		    sym = iupacSym.toStr(m)
## 		except KeyError:
## 		    sym = None
## 		if m in coremannose:
## 		    # Probably N-glycan core Manose
## 		    if c > 3:
## 			bad = True
## 			break
## 		elif c > self.maxlinks.get(sym,2):
## 		    bad = True
## 		    break
## 	    if bad:
## 		continue		
## 	    bad = False
## 	    for m,c in counts.items():
## 		if c > 1:
## 		    bad = True
## 		    break
## 	    if bad:
## 		continue		
## 	    # print counts,counts1
## 	    self.add_instantiation(inst)

    def set_instantiation(self,inst):
        conn = set()
        todo = []
        if self.root():
            todo.append(self.root())
        while len(todo) > 0:
            m = todo.pop(0)
            for l in m.links(False):
                if l.undetermined():
		    if l in inst:
                        l.set_instantiated(True)
			conn.add(l.child())
                todo.insert(0,l.child())
	for ur in self.undetermined_roots():
	    ur.set_connected(ur in conn)
        return

    def instantiations(self):
        if not self._undetermined:
	    yield self
	    return
        plsets = []
        for ur in self.undetermined_roots():
	    if not ur.connected():
                plsets.append(ur.parent_links())
        for inst in combinatorics.product(*plsets,accumulator=combinatorics.set_accumulator):
            self.set_instantiation(inst)
            yield self
        return

    def instantiate(self):
        if not self._undetermined:
            return self
        for g in self.instantiations():
            break
        return self

    def uninstantiate(self):
        if not self._undetermined:
            return self
        self.set_instantiation(set())
        return self
        
    def instantiation_count(self):
        total = 1
        for ur in self.undetermined_roots():
            total *= len(ur.parent_links())
	return total

    def dfsvisit(self,f,m=None,subst=False):
        if m == None:
            self.dfsvisit(f,self.root(),subst)
            for r in self.unconnected_roots():
                self.dfsvisit(f,r,subst)
        else:
            f(m)
            if subst:
                for s in m.substituents():
                    f(s)
            for c in m.children():
                self.dfsvisit(f,c,subst)

    def dfsvisit_post(self,f,m=None,subst=False):
        if m == None:
            self.dfsvisit_post(f,self.root(),subst)
            for r in self.unconnected_roots():
                self.dfsvisit_post(f,r,subst)
        else:
            if subst:
                for s in m.substituents():
                    f(s)
            for c in m.children():
                self.dfsvisit_post(f,c,subst)
            f(m)

    class SubtreeCompositionVisit:
        def __init__(self,sym=None,comp=None):
            
            self.sym = sym
            self.comp = comp
            
        def visit(self,m):

            if self.comp:
                eltcomp = m.composition(self.comp)
                for c in m.children():
                    eltcomp.add(c._elemental_composition)

                m._elemental_composition = eltcomp

            if self.sym:
                symcomp = Composition()
                symcomp[self.sym.toStr(m)] = 1

                for c in m.children():
                    symcomp.add(c._symbol_composition)

                m._symbol_composition = symcomp

    class ElementalCompositionVisit:
        def __init__(self,comp):
            
            self.table = comp
            self.eltcomp = Composition()
            
        def visit(self,m):
            self.eltcomp.add(m.composition(self.table))

    def subtree_composition(self,m,sym_table=None,comp_table=None):
        assert not self.undetermined()
        if m == None:
            m = self.root()
        scv = Glycan.SubtreeCompositionVisit(sym=sym_table,comp=comp_table)
        self.dfsvisit_post(scv.visit,m)

    def elemental_composition(self,comp_table):
        eltcomp = Composition()
        for m in self.all_nodes():
            ec = m.composition(comp_table)
            eltcomp.add(ec)
        return eltcomp

    def byions(self,force=False):
        bions = []
        yions = []
        r = self.root()
        if force or (not hasattr(r,'_symbol_composition') or not hasattr(r,'_elemental_composition')):
            self.subtree_composition(r,sym_table=iupacSym,comp_table=ctable)
        for l in self.all_links():
            # yi,bi = self.split_clone(l)
            c = l.child()
            bions.append((c._symbol_composition,c._elemental_composition,l))
            symcomp = copy.copy(r._symbol_composition)
            symcomp.sub(c._symbol_composition)
            eltcomp = copy.copy(r._elemental_composition)
            eltcomp.sub(c._elemental_composition)
            yions.append((symcomp,eltcomp,l))
        return bions,yions

    def composition(self,force=False):
        r = self.root()
        if force or (not hasattr(r,'_symbol_composition') or not hasattr(r,'_elemental_composition')):
            self.subtree_composition(r,sym_table=iupacSym,comp_table=ctable)
        return r._symbol_composition,r._elemental_composition

    def native_elemental_composition(self):
        return self.elemental_composition(ctable)
    
    def permethlyated_elemental_composition(self):
        return self.elemental_composition(pctable)
    
    def fragments(self,r=None,force=False):
        atroot = False
        if r == None:
            r = self.root()
            atroot = True
            if force or (not hasattr(r,'_symbol_composition') or not hasattr(r,'_elemental_composition')):
                self.subtree_composition(r,sym_table=iupacSym,comp_table=ctable)
        links = r.links()
        nlink = len(links)

        if nlink == 0:
            # self
            fr = (copy.copy(r._symbol_composition),\
                  copy.copy(r._elemental_composition),True,0)
            yield fr
            return
        
        fragstore0 = []
        fragstore1 = []
        for l in links:
            fragstore0.append([])
            fragstore1.append([])
            for fr in self.fragments(l.child()):
                if fr[2]:
                    fragstore0[-1].append(fr)
                else:
                    fragstore1[-1].append(fr)
            fragstore0[-1].append((Composition(),Composition(),True,1))

        for i,prd in enumerate(product(*fragstore0)):
            symcomp = copy.copy(r._symbol_composition)
            eltcomp = copy.copy(r._elemental_composition)
            cl = 0
            for l,fr in zip(links,prd):
                # determine the amount to substract
                symcomp1 = copy.copy(l.child()._symbol_composition)
                symcomp1.sub(fr[0])
                symcomp.sub(symcomp1)
                eltcomp1 = copy.copy(l.child()._elemental_composition)
                eltcomp1.sub(fr[1])
                eltcomp.sub(eltcomp1)
                cl += fr[3]
            fr = (symcomp,eltcomp,True,cl)
            yield fr

        for i in range(nlink):
            for fr in fragstore0[i][:-1]:
                fr1 = (fr[0],fr[1],False,fr[3]+1)
                yield fr1
            for fr in fragstore1[i]:
                fr1 = (fr[0],fr[1],False,fr[3])
                yield fr1

    def subtree_nodes(self,root,subst=False):
        todo = [root]
        seen = set()
        while len(todo) > 0:
            m = todo.pop(0)
            if m not in seen:
                seen.add(m)
                yield m
            if subst:
                for s in m.substituents():
                    if s not in seen:
                        seen.add(s)
                        yield s
            for c in reversed(m.children()):
                todo.insert(0,c)

    def all_nodes(self,subst=False):
        todo = []
        if self.root():
            todo.append(self.root())
	for ur in self.unconnected_roots():
            todo.append(ur)
        for root in todo:
            for m in self.subtree_nodes(root,subst):
                yield m

    def subtree_links(self,root,subst=False,uninstantiated=False):
        for m in self.subtree_nodes(root):
            if subst:
                for sl in m.substituent_links():
                    yield sl
            for l in m.links(instantiated_only=(not uninstantiated)):
                yield l

    def all_links(self,subst=False,uninstantiated=False):
        for m in self.all_nodes():
            if subst:
                for sl in m.substituent_links():
                    yield sl
            for l in m.links(instantiated_only=(not uninstantiated)):
                yield l


    def clone(self):
        self.set_ids()
        if self.root():
            g = Glycan(self.root().deepclone())
        else:
            g = Glycan()
        newurs = set()
        for l in g.all_links(uninstantiated=True):
            if l.undetermined():
                newurs.add(l.child())
        for ur in self.undetermined_roots():
            if len(ur.parent_links()) == 0:
                newurs.add(ur.deepclone())
        g.set_undetermined(newurs)
        return g

    def clone_with_identified_link(self,link):
        assert not self.undetermined()
        r,l = self.root().deepclone(identified_link=link)
        return Glycan(r),l

    def split_clone(self,link):
        g,l = self.clone_with_identified_link(link)
        f = Glycan(l.child())
        l.parent().del_link(l)
        l.child().del_parent_link(l)
        return g,f

    def equals(self,g):

        # Three cases, at least, sigh.

        # 1) no root => compositions, no edges => "trivial" graph
        # isomorphism based on enumeration of potential node matchings

        # 2) tree (not undetermined) => subtree-based

        # 3) rooted DAG (due to undetermined nodes) => graph
        # isomorphism

        # Cases 1 & 3 are implemented by the same code

        # ids in both instances are reset in equal. If returns True,
        # then the ids of each monosaccharide in each glycan will match
        # their counterpart.

        self.set_ids()        
        g.unset_ids()

        if self.has_root() and g.has_root():

            if not self.undetermined() and not g.undetermined():

                # print >>sys.stderr, "Tree comparison"
                g.unset_ids()
                # both are trees, use subtree_equals()
                return self.root().subtree_equals(g.root(),mapids=True)

            else:

                # print >>sys.stderr, "Tree comparison shortcut"
                # both are trees, use subtree_equals()
                if not self.root().subtree_equals(g.root(),mapids=False):
                    return False

        # print >>sys.stderr, "Graph isomorphism comparison"

        g.set_ids()

        nodeset1 = list(self.all_nodes(subst=False))
        nodeset2 = list(g.all_nodes(subst=False))

        if len(nodeset1) != len(nodeset2):
            return False

        linkset1 = set()
        for l in self.all_links(uninstantiated=True):
            linkset1.add((l.parent().id(),l.astuple(),l.child().id()))

        linkset2 = set()
        for l in g.all_links(uninstantiated=True):
            linkset2.add((l.parent().id(),l.astuple(),l.child().id()))

        if len(linkset1) != len(linkset2):
            return False

        # print >>sys.stderr, " ".join(map(lambda i: "%2s"%(i.id(),),nodeset1))
        # print >>sys.stderr, " ".join(map(lambda i: "%2s"%(i.id(),),nodeset2))

        iters = 0
        for ii,jj in iterecmatchings(nodeset1, nodeset2,
                                     self.monosaccharide_match):

            iters += 1
            matching = dict(zip(map(lambda m: m.id(),ii),map(lambda m: m.id(),jj)))
            # print >>sys.stderr, " ".join(map(lambda i: "%2s"%(i.id(),),ii))
            # print >>sys.stderr, " ".join(map(lambda i: "%2s"%(i.id(),),jj))
            linkset3 = set()
            for f,l,t in linkset1:
                linkset3.add((matching[f],l,matching[t]))

            if linkset3 == linkset2:
                for mi,mj in zip(ii,jj):
                    mj.set_id(mi.id())
                # print >>sys.stderr, "%d iterations to find an isomorphism"%(iters,)
                return True

        return False

    @staticmethod
    def monosaccharide_match(a,b):
        # print a
	# print b
        if not a.equals(b):
            return False
        # parent_links_match = False
        # for ii,jj in itermatchings(a.parent_links(),b.parent_links(),
        #                            lambda i,j: i.equals(j) and i.parent().equals(j.parent())):
        #     parent_links_match = True
        #     break
        # if not parent_links_match:
        #     return False
	if len(a.parent_links()) != len(b.parent_links()):
	    return False
	if len(a.links(instantiated_only=True)) != len(b.links(instantiated_only=True)):
	    return False
	if len(a.links(instantiated_only=False)) != len(b.links(instantiated_only=False)):
	    return False
        child_links_match = False
        for ii,jj in itermatchings(a.links(instantiated_only=False),b.links(instantiated_only=False),
                                   lambda i,j: i.equals(j) and i.child().equals(j.child())):
            child_links_match = True
            break
        return child_links_match

    @staticmethod
    def undetroot_equals(a,b,mapids=True):
	if not a.subtree_equals(b,mapids=mapids):
	    return False
	assert(None not in set(l.parent().id() for l in a.parent_links()))
	assert(None not in set(l.parent().id() for l in b.parent_links()))
	uipars = set((l.astuple(),l.parent().id()) for l in a.parent_links())
        ujpars = set((l.astuple(),l.parent().id()) for l in b.parent_links())
        if not (uipars == ujpars):
            return False
	return True

    def str(self,node=None,prefix="",codeprefix="",monofmt=lcSym):
        if node == None:
            node = self.root()
        code = monofmt.toStr(node)
        s = codeprefix + code
        kidlinks = sorted(filter(lambda l: l.instantiated(),node.links()),key=lambda l: Linkage.posstr(l.parent_pos()),reverse=True)
        kids = map(Linkage.child,kidlinks)
        n = len(kids)
        assert n in (0,1,2,3)
        if n == 0:
            return prefix + s
        if n == 1:
            return self.str(kids[0],prefix=prefix,codeprefix=(s + ' - '),monofmt=monofmt)
        if n == 2:
            return self.str(kids[0],prefix + ' '*len(s)+"   ",monofmt=monofmt) + '\n' + \
                   prefix + s + ' + ' + '\n' + \
                   self.str(kids[1],prefix + ' '*len(s)+"   ",monofmt=monofmt)
        if n == 3:
            return self.str(kids[0],prefix + ' '*len(s)+"   ",monofmt=monofmt) + '\n' + \
                   self.str(kids[1],prefix,codeprefix = s + ' + ',monofmt=monofmt) + '\n' + \
                   self.str(kids[2],prefix + ' '*len(s)+"   ",monofmt=monofmt)

    def __str__(self):
        return self.str()

    def dump(self, m=None, level=0, branch='', monofmt=iupacSym):
        if m == None:
            m = self.root()
            
        br = branch + " " + monofmt.toStr(m)
        child_list = []

        for link in m.links():
            child_list.append(link.child())

        if len(child_list) == 0:
            print '    '*level + br
        elif len(child_list) > 1:
            print '    '*level + br
            level += 1
            for c in child_list:
                self.dump(c,level, '', monofmt)
        elif len(child_list) == 1:
            self.dump(child_list[0],level,br,monofmt)

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
    f1 = mf.new('Fuc')
    gc4.add_child(f1,parent_pos=4,child_pos=1,
                  parent_type=Linkage.oxygenPreserved,
                  child_type=Linkage.oxygenLost)
    
    g = Glycan(gc1)
    
    g.dump()

    print g.str(),'\n'
    # bions,yions = g.byions()
    # for bi,yi in zip(bions,yions):
    #    c1,c2,bstr = bi
    #      print "B:",c1,c2,"\n",bstr
    #     c1,c2,ystr = yi
    #    print "Y:",c1,c2,"\n",ystr

    seen = set()
    # for fr in g.fragments():
    for fr in sorted(g.fragments(),key=lambda fr: (fr[1].mass(elmt),fr[3])):
        # if (str(fr[0]),fr[2]) in seen:
        #     continue
        print "%7.2f"%fr[1].mass(elmt),"cl=%d"%fr[3],"Y=%d"%(1*fr[2]),fr[0]
        seen.add((str(fr[0]),fr[2]))

    # for fr in sorted(set(map(lambda fr: str(fr[0]),g.fragments()))):
    #     print fr
                  
                
