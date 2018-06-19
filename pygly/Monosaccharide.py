
import copy
from combinatorics import select

class SuperClass:
    TRI   = 3
    TETRA = 4
    PENT  = 5
    HEX   = 6
    HEPT  = 7
    OCT   = 8
    NON   = 9
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

class Anomer:
    alpha      = 1
    beta       = 2
    uncyclized = 3
    missing    = None

class Monosaccharide:

    # Note that we use None to indicate unset or null values
    # throughout. To unset a value, assign None to it.

    def __init__(self,m=None):

        # Configuration of anomeric carbon 
        self._anomer = None
        
        # Absolute configuration of monosaccharide (D- or L- isomer?)
        # It may be None, representing unset
        self._config = ()

        # Stem 3-letter code [glc, gal, man, rib, gro]
        # It may be None, representing unset
        self._stem = ()

        # Stem and config must have the same number of elements.
        # Could the config be unknown while the Stem is known?
        # Perhaps this should be a list of pairs like mods?

        # Superclass (based on the number of consecutive carbon atoms)
        self._superclass = None

        # Ring closure starting position (THE NUMBER OF THE CARBON) (int)
        self._ring_start = None
            
        # Ring closure ending position (THE NUMBER OF THE CARBON) (int)
        self._ring_end = None

        # Modifier type (Optional, Repeatable)
        self._mods = []

        # Contains a list of links specifying residues that the monosaccharide is linked to
        # Linkage objects in no particular order...
        self._links = []

        # A list of 0 or more links to substituent objects.
        self._substituent_links = []

        # These are not primary data associated with monosaccharides, but are
        # useful for a variety of general processing tasks...
        self._id = None

    def clone(self):
        m = Monosaccharide()
        m._anomer = self._anomer
        m._config = copy.copy(self._config)
        m._stem = copy.copy(self._stem)
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
        # m._id = m._id
        return m

    def noring(self):
        return (self.ring() == (0,0))

    def clear_links(self):
        self._links = []

    def deepclone(self,identified_link=None):
	m = self.clone()
        identified_link_copy = None
	for l in self._links:
	    if identified_link:
	        c,idlc = l.child().deepclone(identified_link=identified_link)
	    else:
		c = l.child().deepclone()
		idlc = None
	    cl = copy.copy(l)
	    cl.set_child(c)
	    cl.set_parent(m)
	    m.add_link(cl)
            if l == identified_link:
                identified_link_copy = cl
            elif idlc != None:
                identified_link_copy = idlc
	if identified_link:
	    return m,identified_link_copy
	return m

    def equals(self,m):
        if self._anomer != m._anomer:
            return False
        if len(self._config) != len(m._config):
            return False
        for a,b in zip(self._config,m._config):
            if a != b:
                return False
        if len(self._stem) != len(m._stem):
            return False
        for a,b in zip(self._stem,m._stem):
            if a != b:
                return False
        if self._superclass != m._superclass:
            return False
        if self._ring_start != m._ring_start:
            return False
        if self._ring_end != m._ring_end:
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

    def anomer(self):
        return self._anomer

    def set_anomer(self,anomer):
        self._anomer = anomer

    def config(self):
        return self._config

    def set_config(self,*config):
        self._config = config

    def stem(self):
        return self._stem

    def set_stem(self,*stem):
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

    def links(self,instantiated_only=True):
        if instantiated_only:
            return filter(lambda l: l.instantiated(),self._links)
        return self._links

    def add_link(self, l):
        self._links.append(l)

    def del_link(self, l):
        self._links.remove(l)

    def id(self):
        return self._id

    def set_id(self,id):
        self._id = id

##     def mass(self,mass_table=None):
##         if not mass_table:
##             mass_table = Monosaccharide.elementMassTable
##         return self._composition.mass(mass_table)

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
	return l

    def has_children(self):
	return (self.first_child() != None)

    def substituents(self):
        return [l.child() for l in self.substituent_links()]

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

        s  = "Monosaccharide:%r\n"%self.id()
        s += "          Anomer = %r\n"%self.anomer()
        s += "          Config = %r\n"%",".join(map(str,self.config()))
        s += "            Stem = %r\n"%",".join(map(str,self.stem()))
        s += "      Superclass = %r\n"%self.superclass()
        s += "            Ring = %r\n"%":".join(map(str,self.ring()))
        if self.has_mods():
            mods = []
            for p,m in self.mods():
                mods.append("%r:%r"%(p,m))
            s += "            Mods = %r\n"%"|".join(mods)
        if self.has_substituents():
            subs = []
            for sl in self.substituent_links():
                subs.append("%s -> %s"%(sl,sl.child()))
            s += "    Substituents = %s\n"%",".join(subs)
        # s += "   isNAcetylated = %s\n"%self.isNAcetylated()
        # s += "     Composition = %s\n"%self.composition()
        # s += "MonoisotopicMass = %f\n"%self.mass()
        if self.has_children():
            s += "      Children = \n        "
            ch = []
            for l in self.links(False):
                if l.instantiated():
                    ch.append("Linkage:%s -> Monosaccharide:%s"%(l,l.child().id()))
                else:
                    ch.append("Linkage:%s ~> Monosaccharide:%s"%(l,l.child().id()))
            s += "\n        ".join(ch) + "\n"
        return s

class Substituent:

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

    def __init__(self,sub):
        self._sub = sub
        self._id = None

    def name(self):
        return self._sub

    def id(self):
        return self._id

    def set_id(self,id):
        self._id = id

    def isNAc(self):
        return (self._sub == Substituent.nAcetyl)

    def has_links(self):
        return False

    def __str__(self):
        if self._id:
            return "%s:%s"%(self._id,self._sub)
        return str(self._sub)

    def equals(self,s):
	return (self._sub == s._sub)

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
        self.set_id(None)
        self.set_parent(None)

    def id(self):
        return self._id

    def set_id(self,id):
        self._id = id

    def child(self):
        return self._child

    def set_child(self, child):
        self._child = child

    def child_type(self):
        return self._child_type

    def set_child_type(self, child_type):
        self._child_type = child_type

    def child_pos(self):
        return self._child_pos

    def set_child_pos(self, child_pos):
        self._child_pos = child_pos

    def child_type2(self):
        return self._child_type2

    def set_child_type2(self, child_type):
        self._child_type2 = child_type

    def child_pos2(self):
        return self._child_pos2

    def set_child_pos2(self, child_pos):
        self._child_pos2 = child_pos

    def parent(self):
        return self._parent

    def set_parent(self, parent):
        self._parent = parent

    def parent_type(self):
        return self._parent_type

    def set_parent_type(self, parent_type):
        self._parent_type = parent_type

    def parent_pos(self):
        return self._parent_pos

    def set_parent_pos(self, parent_pos):
        self._parent_pos = parent_pos

    def parent_type2(self):
        return self._parent_type2

    def set_parent_type2(self, parent_type):
        self._parent_type2 = parent_type

    def parent_pos2(self):
        return self._parent_pos2

    def set_parent_pos2(self, parent_pos):
	# assert(parent_pos == None)
        self._parent_pos2 = parent_pos

    def set_undetermined(self,und):
	self._undetermined = und

    def undetermined(self):
	return self._undetermined

    def set_instantiated(self,inst):
	self._instantiated = inst

    def instantiated(self):
	return self._instantiated

    def compatible(self,a):
	if self._parent_type and a._parent_type and self._parent_type != a._parent_type:
	    return False
	if self._child_type and a._child_type and self._child_type != a._child_type:
	    return False
	if self._child_pos and a._child_pos and self._child_pos != a._child_pos:
	    return False
	ppself = Linkage.parentpos2set(self)
	ppa    = Linkage.parentpos2set(a)
	if len(ppself) > 0 and len(ppa) > 0 and ppself != ppa:
	    return False
	return True

    def compatiblewith(self,a):
	if a._parent_type and self._parent_type != a._parent_type:
	    return False
	if a._child_type and self._child_type != a._child_type:
	    return False
	if a._child_pos and self._child_pos != a._child_pos:
	    return False
	ppself = Linkage.parentpos2set(self)
	ppa    = Linkage.parentpos2set(a)
	if len(ppa) > 0 and ((len(ppself) == 0) or not (ppself <= ppa)):
	    return False
	return True

    @staticmethod
    def parentpos2set(a):
	ppa = set()
	if a._parent_pos:
	    ppa.add(a._parent_pos)
	if a._parent_pos2:
	    ppa.add(a._parent_pos2)
	return ppa

    def equals(self,a):
	if self._parent_type != a._parent_type:
	    return False
	if self._child_type != a._child_type:
	    return False
	if self._child_pos != a._child_pos:
	    return False
	ppself = Linkage.parentpos2set(self)
	ppa    = Linkage.parentpos2set(a)
	if ppself != ppa:
	    return False
	return True

    def __str__(self):
        return "%s(%d+%d)%s"%(self._parent_type if self._parent_type else 'x',
                              self._parent_pos if self._parent_pos else -1,
                              self._child_pos if self._child_pos else -1,
                              self._child_type if self._child_type else 'x')

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

