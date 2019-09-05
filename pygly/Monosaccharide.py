
import copy

from combinatorics import select, itermatchings


class SuperClass:
    TRI   = 3
    TETRA = 4
    PENT  = 5
    HEX   = 6
    HEPT  = 7
    OCT   = 8
    NON   = 9
    DEC   = 10
    missing = None

class Stem:
    gro = 1
    ery = 2
    rib = 3
    ara = 4
    all = 5
    alt = 6
    glc = 7
    man = 8
    tre = 9
    xyl = 10
    lyx = 11
    gul = 12
    ido = 13
    gal = 14
    tal = 15
    thr = 16
    missing = None
    
class Config:
    d = 1
    l = 2
    missing = None

class Mod:
    d       = 1
    keto    = 2
    en      = 3
    a       = 4
    aldi    = 5
    sp2     = 6
    sp      = 7
    geminal = 8
    enx     = 9

class Anomer:
    alpha      = 1
    beta       = 2
    uncyclized = 3
    missing    = None

class Node(object):

    # The Node class represents the functionality common to Monosaccaharides and Substituents...

    def __init__(self):

        # Contains a list of links specifying residues that the monosaccharide is linked to
        # Linkage objects in no particular order...
        self._links = []
        
        # Useful primarily for undetermined roots
        self._parent_links = []

        # Mark undetermined roots that have instantiated links to them...
        self._connected = True

        # These are not primary data associated with monosaccharides, but are
        # useful for a variety of general processing tasks...
        self._id = None

        return 

    def links(self,instantiated_only=True):
        if instantiated_only:
            return filter(lambda l: l.instantiated(),self._links)
        return self._links

    def parent_links(self):
        return self._parent_links

    def add_link(self, l):
        self._links.append(l)

    def del_link(self, l):
        self._links.remove(l)

    def clear_links(self):
        self._links = []

    def add_parent_link(self, l):
        self._parent_links.append(l)

    def del_parent_link(self, l):
        self._parent_links.remove(l)

    def clear_parent_links(self):
        self._parent_links = []

    def children(self):
        return [l.child() for l in self.links() ]

    def first_child(self):
	for l in self.links():
            return l.child()
	return None

    def add_child(self,m,**kw):
        l = Linkage(child=m,**kw)
        self.add_link(l)
        l.set_parent(self)
        m.add_parent_link(l)
	return l

    def has_children(self):
	return (self.first_child() != None)

    def set_connected(self,conn):
	self._connected = conn

    def connected(self):
	return self._connected

    def id(self):
        return self._id

    def set_id(self, id):
        self._id = id

    def unset_id(self):
        self._id = None

    def subtree_equals(self,m,mapids=True):

        if not self.equals(m):
            return False
        
        if mapids:
            m.set_id(self.id())

        for ii,jj in itermatchings(self.links(),m.links(),
                                   lambda i,j: i.equals(j) and i.child().subtree_equals(j.child(),mapids=mapids)):
            return True

        if mapids:
            m.unset_id()

        return False

    def has_substituents(self):
        return False

    def substituents(self):
        return []

    def substituent_links(self):
        return []

    def is_monosaccharide(self):
        return False

class Monosaccharide(Node):

    # Note that we use None to indicate unset or null values
    # throughout. To unset a value, assign None to it.

    def __init__(self):

        super(Monosaccharide,self).__init__()

        # Configuration of anomeric carbon 
        self._anomer = None
        
        # Absolute configuration of monosaccharide (D- or L- isomer?)
        # It may be None, representing unset
        self._config = None

        # Stem 3-letter code [glc, gal, man, rib, gro]
        # It may be None, representing unset
        self._stem = None

        # Superclass (based on the number of consecutive carbon atoms)
        self._superclass = None

        # Ring closure starting position (THE NUMBER OF THE CARBON) (int)
        self._ring_start = None
            
        # Ring closure ending position (THE NUMBER OF THE CARBON) (int)
        self._ring_end = None

        # Modifier type (Optional, Repeatable)
        self._mods = []

        # A list of 0 or more links to substituent objects.
        self._substituent_links = []

        # Short for external descriptor(eg. WURCS) identifier
        # Could be used for keeping track of monosaccharide from descriptor assigned ID
        self._eid = None

    def is_monosaccharide(self):
        return True

    def clone(self):
        m = Monosaccharide()
        m._anomer = self._anomer
        m._config = copy.deepcopy(self._config)
        m._stem = copy.deepcopy(self._stem)
        m._superclass = self._superclass
        m._ring_start = self._ring_start
        m._ring_end = self._ring_end
        m._mods = copy.deepcopy(self._mods)
        # m._composition = copy.copy(self._composition)
        # Will this do it?
        m._substituent_links = copy.deepcopy(self._substituent_links)
        for l in m._substituent_links:
            l.set_parent(m)
        # m._links = copy.deepcopy(self._links)
        m._id = self._id
        m._connected = self._connected
        return m

    def noring(self):
        return (self.ring() == (0,0))

    def deepclone(self,identified_link=None,cache=None):
        if cache == None:
            cache = dict()
	m = self.clone()
        identified_link_copy = None
	for l in self._links:
            assert l.child().id() != None
            if l.child().id() not in cache:
                if identified_link:
                    c,idlc = l.child().deepclone(identified_link=identified_link,cache=cache)
                else:
                    c = l.child().deepclone(cache=cache)
                    idlc = None
                cache[l.child().id()] = c
            else:
                c = cache[l.child().id()]
                idlc = None
	    cl = copy.deepcopy(l)
	    cl.set_child(c)
	    cl.set_parent(m)
	    m.add_link(cl)
            c.add_parent_link(cl)
            if l == identified_link:
                identified_link_copy = cl
            elif idlc != None:
                identified_link_copy = idlc
	if identified_link:
	    return m,identified_link_copy
	return m

    def anysubstmatching(self,m):
        for ii,jj in itermatchings(self.substituent_links(),m.substituent_links(),
                                   lambda i,j: i.equals(j) and i.child().equals(j.child())):
            return True
        return False

    def equals(self,m):
        if not isinstance(m,Monosaccharide):
            return False
        if self._anomer != m._anomer:
            return False
        if self._config != m._config:
            return False
        if self._stem != m._stem:
            return False
        if self._superclass != m._superclass:
            return False
        if self._ring_start != m._ring_start:
            return False
        if self._ring_end != m._ring_end:
            return False
        if self._mods != m._mods:
            return False
        if not self.anysubstmatching(m):
            return False
        return True

    def compatible(self,m):
        if self._anomer and m._anomer and self._anomer != m._anomer:
            return False
        if len(self._config) != len(m._config):
            return False
        for a,b in zip(self._config,m._config):
            if a and b and a != b:
                return False
        if len(self._stem) != len(m._stem):
            return False
        for a,b in zip(self._stem,m._stem):
            if a and b and a != b:
                return False
        if self._superclass and m._superclass and self._superclass != m._superclass:
            return False
        if self._ring_start and m._ring_start and self._ring_start != m._ring_start:
            return False
        if self._ring_end and m._ring_end and self._ring_end != m._ring_end:
            return False
        if len(self._mods) != len(m._mods):
            return False
        for a,b in zip(self._mods,m._mods):
            if a != b:
                return False
        if len(self._substituent_links) != len(m._substituent_links):
            return False
        for a,b in zip(self._substituent_links,m._substituent_links):
            if a.child().name() != b.child().name():
                return False
        return True

    def compatiblewith(self,m,root=False,visibleonly=False):
        if (not visibleonly or not root) and m._anomer and self._anomer != m._anomer:
            return False
	if not visibleonly:
          if len(self._config) != len(m._config):
            return False
          for a,b in zip(self._config,m._config):
            if b and a != b:
                return False
        if len(self._stem) != len(m._stem):
            return False
        for a,b in zip(self._stem,m._stem):
            if b and a != b:
                return False
        if m._superclass and self._superclass != m._superclass:
            return False
        if not visibleonly and m._ring_start and self._ring_start != m._ring_start:
            return False
        if not visibleonly and m._ring_end and self._ring_end != m._ring_end:
            return False
        if len(self._mods) != len(m._mods):
            return False
        for a,b in zip(self._mods,m._mods):
            if a != b:
                return False
        if len(self._substituent_links) != len(m._substituent_links):
            return False
        anymatch = False
	for slinks in select(m._substituent_links,len(self._substituent_links)):
	  badmatch = False
          for a,b in zip(self._substituent_links,slinks):
	    if not a.equals(b) or not a.child().equals(b.child()):
                badmatch = True
	  if not badmatch:
	    anymatch = True
	    break
	if not anymatch:
	  return False
        return True

    def fully_determined(self):
        if self._anomer == Anomer.missing:
            return False
        if self._config == Config.missing:
            return False
        if self._stem == Stem.missing:
            return False
        if self._superclass == SuperClass.missing:
            return False
        if self._ring_start == None:
            return False
        if self._ring_end == None:
            return False
        for m in self._mods:
            if m[0] == None:
                return False
        return True

    def anomer(self):
        return self._anomer

    def set_anomer(self,anomer):
        self._anomer = anomer

    def config(self):
        return self._config

    def set_config(self,*config):
	if set([config]) == set([Config.missing]):
	    self._config = None
        else:
            self._config = config

    def stem(self):
        return self._stem

    def set_stem(self,*stem):
	if set([stem]) == set([Stem.missing]):
	    self._stem = None
	else:
            self._stem = stem

    def superclass(self):
        return self._superclass

    def set_superclass(self,cls):
        self._superclass = cls

    def ring_start(self):
        return self._ring_start

    def set_ring_start(self,start):
        self._ring_start = start

    def ring_end(self):
        return self._ring_end

    def set_ring_end(self,end):
        self._ring_end = end

    def ring(self):
        return (self.ring_start(),self.ring_end())

    def mods(self):
        return self._mods

    def add_mod(self,pos,mod):
        if isinstance(pos,str):
            pos = tuple(sorted(map(int,pos.split(','))))
        else:
            pos = (int(pos),)
        self._mods.append((pos,mod))
        self._mods = sorted(self._mods)

    def remove_mod(self,mod,pos=None):
        for i in range(len(self._mods)-1,-1,-1):
            if (self._mods[i][1]) == mod and (pos == None or self._mods[i][0] == pos):
                del self._mods[i]

    def count_mod(self,mod=None):
        count = 0
        for m in self._mods:
            if m[1] == mod or mod == None:
                count += 1
        return count

    def clear_mods(self):
        self._mods = []

    def has_mods(self):
        return len(self._mods) > 0

    def composition(self,comp_table):
        c = comp_table.new()
        if self.superclass() != None:
            c.add(comp_table[('SuperClass',self.superclass())])
        for pos,mod in self.mods():
            c.add(comp_table[('Mod',mod)])
        for sub in self.substituents():
            c.add(comp_table[('Substituent',sub.name())])
	return c

##     def set_composition(self,**kw):
##         self._composition.update(kw)

    def substituent_links(self,instantiated_only=True):
        if instantiated_only:
            return filter(lambda l: l.instantiated(),self._substituent_links)
        return self._substituent_links

    def add_substituent_link(self, l):
        self._substituent_links.append(l)

    def set_external_descriptor_id(self, eid):
        self._eid = eid

    def external_descriptor_id(self):
        return self._eid

##     def mass(self,mass_table=None):
##         if not mass_table:
##             mass_table = Monosaccharide.elementMassTable
##         return self._composition.mass(mass_table)

    def substituents(self):
        return set([l.child() for l in self.substituent_links()])

    def add_substituent(self,sub,**kw):
        if isinstance(sub,Substituent):
            l = SubLinkage(child=sub,**kw)
        else:
            l = SubLinkage(child=Substituent(sub),**kw)
        self.add_substituent_link(l)
        l.set_parent(self)
	return l

    def has_substituents(self):
        return len(self._substituent_links) > 0

    def is_nacetylated(self):
	sl = self.substituents()
	if len(sl) != 1:
	    return False
	return sl[0].isNAc()

    def has_non_nacetyl_substituents(self):
        for s in self.substituents():
            if not s.isNAc():
                return True
        return False

    def __str__(self):

        s  = "Monosaccharide:%s\n"%self.id()
        s += "          Anomer = %s\n"%constantString(Anomer,self.anomer())
        s += "          Config = %s\n"%constantStrings(Config,self.config())
        s += "            Stem = %s\n"%constantStrings(Stem,self.stem())
        s += "      Superclass = %s\n"%constantString(SuperClass,self.superclass())
        s += "            Ring = %s\n"%":".join(map(str,self.ring()))
        if self.has_mods():
            mods = []
            for p,m in self.mods():
                mods.append("%s:%s"%(",".join(map(str,p)),constantString(Mod,m)))
            s += "            Mods = %s\n"%"|".join(mods)
        if self.has_substituents():
            subs = []
            for sl in self.substituent_links():
                subs.append("%s -> %s"%(sl,sl.child()))
            s += "    Substituents = %s\n"%"\n                   ".join(subs)
        # s += "   isNAcetylated = %s\n"%self.isNAcetylated()
        # s += "     Composition = %s\n"%self.composition()
        # s += "MonoisotopicMass = %f\n"%self.mass()
        if self.has_children():
            s += "        Children = "
            ch = []
            for l in self.links(False):
                if l.instantiated():
                    ch.append("%s -> Monosaccharide:%s"%(l,l.child().id()))
                else:
                    ch.append("%s ~> Monosaccharide:%s"%(l,l.child().id()))
            s += "\n                   ".join(ch) + "\n"
        if len(self.parent_links()) > 0:
            s += "         Parents = "
            ch = []
            for l in self.parent_links():
                if l.instantiated():
                    ch.append("Monosaccharide:%s -> %s"%(l.parent().id(),l))
                else:
                    ch.append("Monosaccharide:%s ~> %s"%(l.parent().id(),l))
            s += "\n                   ".join(ch) + "\n"

        return s

class Substituent(Node):

    #Substituents
    nAcetyl             = 1
    sulfate             = 2
    amino               = 3
    phosphoethanolamine = 4
    nglycolyl           = 5
    acetyl              = 6
    nsulfate            = 7
    anhydro             = 8
    phosphate           = 9
    methyl              = 10
    nmethyl             = 11
    nsuccinate          = 12
    pyruvate            = 13
    diphosphoethanolamine = 14
    pyrophosphate       = 15
    formyl              = 16
    thio                = 17
    fluoro		= 18
    bromo		= 19
    glycolyl		= 20
    ndimethyl		= 21
    chloro		= 22
    ethyl		= 23
    succinate		= 24
    nformyl		= 25
    spyruvate		= 26
    rpyruvate		= 27
    namidino		= 28
    iodo		= 29
    ethanolamine	= 30 
    slactate		= 31
    rlactate		= 32
    xlactate		= 33
    hydroxymethyl	= 34
    triphosphate	= 35
    lactone		= 36
    rcarboxyethyl	= 37
    scarboxyethyl	= 38
    phosphocholine	= 39
    methyl_oxygen_lost  = 40
    sulfate_oxygen_lost = 41
    amino_oxygen_preserved = 42
    phosphate_oxygen_lost = 43
    acetyl_oxygen_lost = 44

    def __init__(self,sub):

        super(Substituent,self).__init__()

        self._sub = sub

    def name(self):
        return self._sub

    def isNAc(self):
        return (self._sub == Substituent.nAcetyl)

    def composition(self,comp_table):
        c = comp_table.new()
        c.add(comp_table[('Substituent',self.name())])
        return c

    def __str__(self):
        if self._id:
            return "%s:%s"%(self._id,constantString(Substituent,self._sub))
        return constantString(Substituent,self._sub)

    def equals(self,s):
        if not isinstance(s,Substituent):
            return False
	return (self._sub == s._sub)

    def fully_determined(self):
        return True

class Linkage:

    # Atom Replacement constants (link types)
    oxygenPreserved     = 1
    oxygenLost          = 2
    hydrogenLost        = 3
    nitrogenAdded       = 4
    missing             = None

    # Note that we use None to indicate unset values throughout. To
    # unset a value, assign None to it.

    def __init__(self, child,
                 parent_type=None, parent_pos=None,
                 child_type=None, child_pos=None,
                 parent_type2=None, parent_pos2=None,
                 child_type2=None, child_pos2=None,
		 undetermined=False):

        # the monosaccharide child in the linkage
        self.set_child(child)
        
        # parent linkage type
        self.set_parent_type(parent_type)

        # position of linkage attachment to parent
        self.set_parent_pos(parent_pos)
        
        # child linkage type
        self.set_child_type(child_type)
        
        # position of linkage attachment to child
        self.set_child_pos(child_pos)

        # 2nd parent linkage type
        self.set_parent_type2(parent_type2)

        # position of 2nd linkage attachment to parent
        self.set_parent_pos2(parent_pos2)
        
        # 2nd child linkage type
        self.set_child_type2(child_type2)
        
        # position of 2nd linkage attachment to child
        self.set_child_pos2(child_pos2)

	# Undetermined status of the link 
	self.set_undetermined(undetermined)
	if undetermined:
	    self.set_instantiated(False)
	else:
	    self.set_instantiated(True)

        # These are not primary data associated with linkages, but are
        # useful for a variety of general processing tasks...
        self.unset_id()
        self.set_parent(None)

    def id(self):
        return self._id

    def reverse(self):
	assert self.parent()
	l = Linkage(child=self.parent(),
                    parent_pos=self.child_pos(),
                    parent_type=self.child_type(),
                    child_pos=self.parent_pos(),
                    child_type=self.parent_type(),
                    undetermined=self.undetermined())
	l.set_parent(self.child())
	self.child().del_parent_link(self)
	self.child().add_link(l)
	self.parent().del_link(self)
	self.parent().add_parent_link(l)
	return l

    def clone(self):
        l = Linkage(child=self.child(),
                    parent_pos=self.parent_pos(),
                    parent_type=self.parent_type(),
                    child_pos=self.child_pos(),
                    child_type=self.child_type(),
                    undetermined=self.undetermined())
        l.set_id(self.id())
        l.set_parent(self.parent())
        l.set_instantiated(self.instantiated())
        return l

    def set_id(self,id):
        self._id = id

    def unset_id(self):
        self._id = None

    def child(self):
        return self._child

    def set_child(self, child):
        self._child = child

    def child_type(self):
        return self._child_type

    def set_child_type(self, child_type):
        # instance setting 
        if child_type in (1,2,3,4):
            self._child_type = set([child_type])
        else:
            # None or set of 1,2,3,4 hopefully
            self._child_type = child_type

    def set_child_type2(self, child_type):
        if child_type == None:
            return
        assert child_type in (1,2,3,4)
        assert self._child_type
        self._child_type.add(child_type)

    def child_pos(self):
        return self._child_pos

    def set_child_pos(self, child_pos):
        try:
            child_pos = int(child_pos)
            if child_pos == -1:
                self._child_pos = None
            else:
                self._child_pos = set([child_pos])
        except (ValueError,TypeError):
            if child_pos == None:
                self._child_pos = None
            else:
                self._child_pos = set(child_pos)

    def set_child_pos2(self, child_pos):
        if child_pos == None:
            return
        assert isinstance(child_pos,int)
        assert self._child_pos
        self._child_pos.add(child_pos)

    def parent(self):
        return self._parent

    def set_parent(self, parent):
        self._parent = parent

    def parent_type(self):
        return self._parent_type

    def set_parent_type(self, parent_type):
        if parent_type in (1,2,3,4):
            # instance setting 
            self._parent_type = set([parent_type])
        else:
            # None or set of 1,2,3,4 hopefully
            self._parent_type = parent_type

    def set_parent_type2(self, parent_type):
        if parent_type == None:
            return
        assert parent_type in (1,2,3,4)
        assert self._parent_type
        self._parent_type.add(parent_type)

    def parent_pos(self):
        return self._parent_pos

    def set_parent_pos(self, parent_pos):
        try:
            parent_pos = int(parent_pos)
            if parent_pos == -1:
                self._parent_pos = None
            else:
                self._parent_pos = set([parent_pos])
        except (ValueError,TypeError):
            if parent_pos == None:
                self._parent_pos = None
            else:
                self._parent_pos = set(parent_pos)

    def set_parent_pos2(self, parent_pos):
        if parent_pos == None:
            return
        assert isinstance(parent_pos,int)
        assert self._parent_pos
        self._parent_pos.add(parent_pos)

    def set_undetermined(self,und):
	self._undetermined = und

    def undetermined(self):
	return self._undetermined

    def set_instantiated(self,inst):
	self._instantiated = inst

    def instantiated(self):
	return self._instantiated

##     def compatible(self,a):
## 	if self._parent_type and a._parent_type and self._parent_type != a._parent_type:
## 	    return False
## 	if self._child_type and a._child_type and self._child_type != a._child_type:
## 	    return False
## 	if self._child_pos and a._child_pos and self._child_pos != a._child_pos:
## 	    return False
## 	ppself = self.parentpos2set()
## 	ppa    = a.parentpos2set()
## 	if len(ppself) > 0 and len(ppa) > 0 and ppself != ppa:
## 	    return False
## 	return True

##     def compatiblewith(self,a):
## 	if a._parent_type and self._parent_type != a._parent_type:
## 	    return False
## 	if a._child_type and self._child_type != a._child_type:
## 	    return False
## 	if a._child_pos and self._child_pos != a._child_pos:
## 	    return False
## 	ppself = self.parentpos2set()
## 	ppa    = a.parentpos2set()
## 	if len(ppa) > 0 and ((len(ppself) == 0) or not (ppself <= ppa)):
## 	    return False
## 	return True

    def equals(self,a):
	# print "---"
	# print repr(self._parent_type), repr(self._child_type), repr(self._child_pos), repr(self._parent_pos)
	# print repr(a._parent_type), repr(a._child_type), repr(a._child_pos), repr(a._parent_pos)
	# print "---"
        # if self._undetermined != a._undetermined:
        #     return False
        if self._instantiated != a._instantiated:
            return False
	if self._parent_type != a._parent_type:
	    return False
	if self._child_type != a._child_type:
	    return False
	if self._child_pos != a._child_pos:
	    return False
        if self._parent_pos != a._parent_pos:
	    return False
	return True

    def fully_determined(self):
        if not self._instantiated:
            return False
        if self._parent_type == None:
	    return False
        if self._child_type == None:
            return False
        if self._parent_pos == None:
            return False
        if self._child_pos == None:
            return False
        return True

    @staticmethod
    def valtuple(val):
        return (tuple(sorted(val)) if val else None)

    @staticmethod
    def typestr(val,delim="|"):
        return (delim.join(map(lambda v: constantString(Linkage,v),Linkage.valtuple(val))) if val else "missing")
                            
    @staticmethod
    def posstr(val,delim="|"):
        return (delim.join(map(str,Linkage.valtuple(val))) if val else -1)
                            
    def astuple(self):
        return (self._instantiated,
                self.valtuple(self._parent_type),
                self.valtuple(self._parent_pos),
                self.valtuple(self._child_pos),
                self.valtuple(self._child_type))

    def __str__(self):
        return "%s (%s+%s) %s"%(self.typestr(self._parent_type),
                                self.posstr(self._parent_pos),
                                self.posstr(self._child_pos),
                                self.typestr(self._child_type))

# Should we specialize substituent linkages?
class SubLinkage(Linkage):
    pass

# This is here to put peptide mass and elemental composition on the same
# footing as glycans....

class AminoAcid:
    A = 1
    C = 2
    D = 3
    E = 4
    F = 5
    G = 6
    H = 7
    I = 8
    K = 9
    L = 10
    M = 11
    N = 12
    P = 13
    Q = 14
    R = 15
    S = 16
    T = 17
    V = 18
    W = 19
    Y = 20

class PeptideTerminal:

    nTerm = 1
    cTerm = 2

class PeptidePTM:

    Carbamidomethyl = 1
    Oxidation = 2
    Phosphorylation = 3
    Methylation = 4
    Dimethylation = 5
    Deamidation = 6
    TMTSixplex = 7

def constantLookup(s):
    cls = s.split('.')[0]
    return (cls,eval(s))

def constantString(cls,i):
    for k in dir(cls):
        if not k.startswith('__') and getattr(cls,k) == i:
            return k
    return i

def constantStrings(cls,i,delim=","):
    try:
        return delim.join(map(lambda ii: constantString(cls,ii),i))
    except TypeError:
        pass
    return constantString(cls,i)

if __name__ == '__main__':

    bdGal = Monosaccharide()

    bdGal.set_id("bdGal")
    bdGal.set_anomer(Anomer.beta)
    bdGal.set_config(Config.d)
    bdGal.set_stem(Stem.gal)
    bdGal.set_superclass(SuperClass.HEX)
    bdGal.set_ring_start(1)
    bdGal.set_ring_end(5)

    print bdGal

    bdGlc = bdGal.clone()
    bdGlc.set_id("bdGlc")
    bdGlc.set_stem(Stem.glc)

    print bdGlc
    print bdGal
    
    #Make GlcNAc...
    GlcNAc = bdGlc.clone()
    GlcNAc.set_id("GlcNAc")
    GlcNAc.add_substituent(Substituent.nAcetyl,
                           parent_pos=2,
                           parent_type=Linkage.oxygenLost,
                           child_type=Linkage.nitrogenAdded)

    print GlcNAc
    print bdGal
    print bdGlc

    #Neu5Ac test
    Neu5Ac = Monosaccharide()

    Neu5Ac.set_id("Neu5Ac")
    Neu5Ac.set_anomer(Anomer.alpha)
    Neu5Ac.set_config(Config.d,Config.d)
    Neu5Ac.set_stem(Stem.gro,Stem.gal)
    Neu5Ac.set_superclass(SuperClass.NON)
    Neu5Ac.set_ring_start(2)
    Neu5Ac.set_ring_end(6)
    Neu5Ac.add_mod(1,Mod.a)
    Neu5Ac.add_mod(2,Mod.keto)
    Neu5Ac.add_mod(3,Mod.d)
    Neu5Ac.add_substituent(Substituent.nAcetyl,
                           parent_pos=5,
                           parent_type=Linkage.oxygenLost,
                           child_pos=1,
                           child_type=Linkage.nitrogenAdded)

    #print test case
    print Neu5Ac

    m = Neu5Ac.clone()
    m.set_id("Clone of Neu5Ac")

    print m

    bdGal.add_child(Neu5Ac,
                    parent_type=Linkage.oxygenPreserved,
                    parent_pos=2,
                    child_pos=2,
                    child_type=Linkage.oxygenLost)

    print bdGal

