

import copy

from . combinatorics import select, itermatchings
try:
    from past.builtins import basestring    
except ImportError:
    pass

class SuperClass:
    TRI   = 3
    TETRA = 4
    PENT  = 5
    HEX   = 6
    HEPT  = 7
    OCT   = 8
    NON   = 9
    DEC   = 10
    SUG   = 20
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

        # Contains a list of "normal" links specifying residues that the monosaccharide is linked to
        # Linkage objects in no particular order...
        self._links = []

        # Special links (uninstantiated or repeat links), split out for efficiency
        self._special_links = []
        
        # Useful primarily for undetermined roots
        self._parent_links = []

        # Mark undetermined roots that have instantiated links to them...
        self._connected = True

        # These are not primary data associated with monosaccharides, but are
        # useful for a variety of general processing tasks...
        self._id = None

        self._external_descriptor = None
        self._eid = None

        return 

    def is_default_link_type(self,l):
        # NON_REPEAT + REPEAT_EXIT, INSTANTIATED
        return l.matching_link_type(None,-Linkage.REPEAT_BRIDGE, Linkage.INSTANTIATED)

    def links(self, default=True, **kw):
        if default:
            for l in self._links:
                yield l
        else:
            for l in self._links + self._special_links:
                if l.matching_link_type(**kw):
                    yield l

    def links_without_repeat(self):
        # NON_REPEAT, INSTANTIATED
        return self.links(default=False,repeat=Linkage.NON_REPEAT,inst=Linkage.INSTANTIATED)

    def links_with_repeat_exit(self):
        # NON_REPEAT + REPEAT_EXIT, INSTANTIATED; same as default
        return self.links(default=False,repeat=-Linkage.REPEAT_BRIDGE,inst=Linkage.INSTANTIATED)

    def links_has_repeat_bridge(self):
        # REPEAT_BRIDGE, INSTANTIATED only
        return any(self.links(default=False,repeat=Linkage.REPEAT_BRIDGE,inst=Linkage.INSTANTIATED))

    def links_repeat_bridge_only(self):
        # REPEAT_BRIDGE, INSTANTIATED only
        return self.links(default=False,repeat=Linkage.REPEAT_BRIDGE,inst=Linkage.INSTANTIATED)

    def links_with_repeat_bridge(self):
        # NON_REPEAT + REPEAT_BRIDGE, INSTANTIATED, similar to default
        return self.links(default=False,repeat=-Linkage.REPEAT_EXIT,inst=Linkage.INSTANTIATED)

    def links_with_uninstantiated(self):
        # NON_REPEAT + REPEAT_EXIT
        return self.links(default=False,repeat=-Linkage.REPEAT_BRIDGE)

    def link_count(self, **kw):
        return sum(1 for _ in self.links(**kw))

    def has_links(self, **kw):
        for l in self.links(**kw):
            return True
        return False

    def any_link(self, **kw):
        for l in self.links(**kw):
            return l
        return None

    def add_link(self, l):
        if self.is_default_link_type(l):
            self._links.append(l)
        else:
            self._special_links.append(l)

    def del_link(self, l):
        if self.is_default_link_type(l):
            self._links.remove(l)
        else:
            self._special_links.remove(l)

    def clear_links(self):
        self._links = []
        self._special_links = []

    def parent_links(self, default=True, **kw):
        if default:
            for l in self._parent_links:
                if l.matching_link_type(None,-Linkage.REPEAT_BRIDGE,None):
                    yield l
        else:
            for l in self._parent_links:
                if l.matching_link_type(**kw):
                    yield l

    def parent_link_count(self,**kw):
        return sum(1 for _ in self.parent_links(**kw))

    def has_parent_links(self,**kw):
        for l in self.parent_links(**kw):
            return True
        return False

    def any_parent_link(self,**kw):
        for l in self.parent_links(**kw):
            return l
        return None

    def add_parent_link(self, l):
        self._parent_links.append(l)

    def del_parent_link(self, l):
        self._parent_links.remove(l)

    def clear_parent_links(self):
        self._parent_links = []

    def parents(self):
        return [ l.parent() for l in self.parent_links() ]

    def children(self):
        return [ l.child() for l in self.links() ]

    def any_child(self):
        for l in self.links():
            return l.child()
        return None

    def add_child(self,m,**kw):
        l = Linkage(child=m,parent=self,**kw)
        return self.add_linkage(l)

    def add_child_with_special_linkage(self, m, l):
        l.set_child(m)
        l.set_parent(self)
        return self.add_linkage(l)

    def add_linkage(self, l):
        self.add_link(l)
        l.child().add_parent_link(l)
        return l

    def has_children(self):
        return (self.any_link() != None)

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

    def is_substituent(self):
        return False

    def external_descriptor(self):
        return self._external_descriptor

    def set_external_descriptor(self, s):
        self._external_descriptor = s

    def set_external_descriptor_id(self, eid):
        self._eid = eid

    def external_descriptor_id(self):
        return self._eid

    def is_repeat_start(self):
        for l in self.parent_links(default=False,repeat=Linkage.REPEAT_BRIDGE,inst=Linkage.INSTANTIATED):
            return True
        return False

    def is_repeat_end(self):
        return self.links_has_repeat_bridge()

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

        # If any substituent has a link to a monosaccharide, it must be indicated by this boolean.
        # Updated by Substituent.add_link
        self._has_subst_out_link = False

    def is_monosaccharide(self):
        return True

    def clone(self):
        # Clones the monosaccharide (and its substituents and substituent links, if necessary)
        # Note that _links etc. are *not* instantiated here (left empty). 
        m = Monosaccharide()
        m._anomer = self._anomer
        m._config = copy.deepcopy(self._config)
        m._stem = copy.deepcopy(self._stem)
        m._superclass = self._superclass
        m._ring_start = self._ring_start
        m._ring_end = self._ring_end
        m._mods = copy.deepcopy(self._mods)
        # Note! This is only valid because substituents have only one parent. 
        seen_subst = dict()
        for sl in self.substituent_links():
            slc = sl.clone()
            if sl.child() in seen_subst:
                ch = seen_subst[sl.child()]
            else:
                ch = sl.child().clone()
                seen_subst[sl.child()] = ch
            slc.set_parent(m)
            slc.set_child(ch)
            m.add_substituent_link(slc)
            ch.add_parent_link(slc)
        # m._composition = copy.copy(self._composition)
        m._id = self._id
        m._connected = self._connected
        m._external_descriptor = self._external_descriptor
        m._eid = self._eid
        return m

    def noring(self):
        return (self.ring() == (0,0))

    def deepclone(self, identified_link=None, cache=None):

        if cache == None:
            cache = dict()

        m = self.clone()
        identified_link_copy = None
        for l in self.links(default=False):
            assert l.child().id() != None
            assert l.child().is_monosaccharide()
            if l.child().id() not in cache:
                # child of repeat bridge link should have already been made
                # all others may be new
                assert not l.is_repeat_bridge()
                if identified_link:
                    c,idlc = l.child().deepclone(identified_link=identified_link,cache=cache)
                else:
                    c = l.child().deepclone(cache=cache)
                    idlc=None
                cache[l.child().id()] = (c,idlc)
            else:
                c,idlc = cache[l.child().id()]

            cl = l.clone()
            cl.set_child(c)
            if not cl.is_subst_out_link():
                cl.set_parent(m)
                m.add_linkage(cl)
            else:
                assert(l.parent().id())
                subst_parent_id = l.parent().id()
                for subst in m.substituents():
                    assert(subst.id())
                    if subst.id() == subst_parent_id:
                        cloned_subst_parent = subst
                        break
                cl.set_parent(cloned_subst_parent)
                cloned_subst_parent.add_linkage(cl)
                
            if l == identified_link:
                identified_link_copy = cl
            elif idlc != None:
                identified_link_copy = idlc

        if identified_link:
            return m, identified_link_copy

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
        if any(map(lambda c: c == Config.missing,self._config)):
            return False
        if self._stem == Stem.missing:
            return False
        if any(map(lambda s: s == Stem.missing,self._stem)):
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

    def root_partially_determined(self):
        # Don't test the ring or anomer configuration
        if self._config == Config.missing:
            return False
        if self._stem == Stem.missing:
            return False
        if self._superclass == SuperClass.missing:
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
        if set(config) == set([Config.missing]):
            self._config = None
        else:
            self._config = config

    def stem(self):
        return self._stem

    def set_stem(self,*stem):
        if set(stem) == set([Stem.missing]):
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
        if isinstance(pos,basestring):
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

    def substituent_links(self):
        for l in self._substituent_links:
            yield l

    def add_substituent_link(self, l):
        self._substituent_links.append(l)

    def remove_substituent_link(self, l):
        self._substituent_links.remove(l)

##     def mass(self,mass_table=None):
##         if not mass_table:
##             mass_table = Monosaccharide.elementMassTable
##         return self._composition.mass(mass_table)

    def substituents(self):
        seen_subst = set()
        for l in self.substituent_links():
            if l.child() not in seen_subst:
                yield l.child()
                seen_subst.add(l.child())

    def substituent_count(self):
        return sum(1 for _ in self.substituents())

    def add_substituent(self,sub,**kw):
        if isinstance(sub,Substituent):
            l = SubLinkage(child=sub,parent=self,**kw)
        else:
            l = SubLinkage(child=Substituent(sub),parent=self,**kw)
        self.add_substituent_link(l)
        l.child().add_parent_link(l)
        return l

    def has_substituents(self):
        return len(self._substituent_links) > 0

    def is_nacetylated(self):
        firstsub = None
        for i,s in enumerate(self.substituents()):
            if i > 0:
                return False
            firstsub = s
        return firstsub != None and firstsub.isNAc()

    def has_non_nacetyl_substituents(self):
        for s in self.substituents():
            if not s.isNAc():
                return True
        return False

    # Monosaccharide links can be: basic, uninstantiated, repeat_exit, repeat_bridge, subst_out
    # By default: basic and subst_out are generated...
    def links(self, **kw):
        for l in super(Monosaccharide,self).links(**kw):
            yield l
        if self._has_subst_out_link:
            for subst in self.substituents():
                for l in subst.links(**kw):
                    yield l

    def __str__(self):

        s  = "Monosaccharide:%s (%r)\n"%(self.id(),self)
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
            for l in self.links(default=False): #implies all links
                l_linked_by_sub = not l.parent().is_monosaccharide()
                l_linkage_symbol = "~"
                if l.instantiated():
                    l_linkage_symbol = "-"
                if l.is_repeat_bridge():
                    l_linkage_symbol = "-[R]-"

                if l_linked_by_sub:
                    l_upper = l.parent().any_parent_link()
                    l_lower = l
                    tmp_s = ""
                    if l.is_repeat_bridge():
                        tmp_s = "[R]"
                    ch.append("%s -( %s )%s-> %s Monosaccharide:%s" % (l_upper, l.parent(), tmp_s, l, l.child().id()))
                else:
                    l_linkage_symbol += ">"
                    ch.append("%s %s Monosaccharide:%s" % (l, l_linkage_symbol, l.child().id()))

                #if l.instantiated():
                #    ch.append("%s -> Monosaccharide:%s"%(l,l.child().id()))
                #else:
                #    ch.append("%s ~> Monosaccharide:%s"%(l,l.child().id()))
            s += "\n                   ".join(ch) + "\n"
        if self.has_parent_links():
            s += "         Parents = "
            ch = []
            for l in self.parent_links():
                l_parent_type = "Substituent"
                if l.parent().is_monosaccharide():
                    l_parent_type = "Monosaccharide"
                else:
                    l_parent_type = "Monosaccharide:%s - *%s(sub)*" % (l.parent().any_parent_link().parent().id(), str(l.parent()) )

                l_linkage_symbol = "~>"
                if l.instantiated():
                    l_linkage_symbol = "->"
                if l.is_repeat_bridge():
                    l_linkage_symbol = "-[R]->"
                ch.append("%s:%s %s %s" % (l_parent_type, l.parent().id(), l_linkage_symbol, l))


                #if l.instantiated():
                #    ch.append("Monosaccharide:%s -> %s"%(l.parent().id(),l))
                #else:
                #    ch.append("Monosaccharide:%s ~> %s"%(l.parent().id(),l))
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
    fluoro              = 18
    bromo               = 19
    glycolyl            = 20
    ndimethyl           = 21
    chloro              = 22
    ethyl               = 23
    succinate           = 24
    nformyl             = 25
    spyruvate           = 26
    rpyruvate           = 27
    namidino            = 28
    iodo                = 29
    ethanolamine        = 30 
    slactate            = 31
    rlactate            = 32
    xlactate            = 33
    hydroxymethyl       = 34
    triphosphate        = 35
    lactone             = 36
    rcarboxyethyl       = 37
    scarboxyethyl       = 38
    phosphocholine      = 39
    methyl_oxygen_lost  = 40
    sulfate_oxygen_lost = 41
    amino_oxygen_preserved = 42
    phosphate_oxygen_lost = 43
    acetyl_oxygen_lost = 44
    phosphate_bridged = 45

    def __init__(self,sub):

        super(Substituent,self).__init__()
        self._id = None
        self._sub = sub

    def clone(self):
        s = Substituent(self.name())
        s.set_id(self.id())
        s.set_connected(self.connected())
        s.set_external_descriptor(self.external_descriptor())
        s.set_external_descriptor_id(self.external_descriptor_id())
        return s

    def name(self):
        return self._sub

    def is_substituent(self):
        return True

    def isNAc(self):
        return (self._sub == Substituent.nAcetyl)

    def composition(self,comp_table):
        c = comp_table.new()
        c.add(comp_table[('Substituent',self.name())])
        return c

    def add_link(self, l):
        # if we ever add a link, we need to update the monosaccharide too
        super(Substituent,self).add_link(l)
        assert(l.matching_link_type(fromto=Linkage.SUBST_TO_MONO))
        for p in self.parents():
            # should be exactly one parent, the monosaccharide this
            # substituent is considered part of
            p._has_subst_out_link = True

    def __str__(self):
        if self._id:
            return "%s:%s (%r)"%(self._id,constantString(Substituent,self._sub),self)
        return "%s (%r)"%(constantString(Substituent,self._sub),self)

    def equals(self,s):
        if not isinstance(s,Substituent):
            return False
        return (self._sub == s._sub)

    def fully_determined(self):
        return True

class Linkage(object):

    # Atom Replacement constants (link types)
    oxygenPreserved     = 1
    oxygenLost          = 2
    hydrogenLost        = 3
    nitrogenAdded       = 4
    missing             = None

    # Note that we use None to indicate unset values throughout. To
    # unset a value, assign None to it.

    # Link types - 3/4 values 
    MONO_TO_MONO = 101
    MONO_TO_SUBST = 102
    SUBST_TO_MONO = 103
    # SUBST_TO_SUBST = 4

    # Repeat status - 3 values
    NON_REPEAT = 201
    REPEAT_BRIDGE = 202
    REPEAT_EXIT = 203

    # Instantiated - 2 values
    INSTANTIATED = 301
    UNINSTANTIATED = 302
    
    # Most links are of this type # Static data-member
    _link_type = (MONO_TO_MONO,NON_REPEAT,INSTANTIATED)

    def __init__(self, child=None, parent=None,
                 parent_type=None, parent_pos=None,
                 child_type=None, child_pos=None,
                 parent_type2=None, parent_pos2=None,
                 child_type2=None, child_pos2=None):

        # the monosaccharide child in the linkage
        self.set_child(child)

        self.set_parent(parent)
        
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

        # These are not primary data associated with linkages, but are
        # useful for a variety of general processing tasks...
        self.unset_id()

    def id(self):
        return self._id

    def link_type(self):
        return self._link_type

    def is_link_type(self,link_type):
        return self._link_type == link_type

    def matching_link_type(self,fromto=None,repeat=None,inst=None):
        if fromto != None:
            if fromto > 0:
                if self._link_type[0] != fromto:
                    return False
            else:
                if self._link_type[0] == -fromto:
                    return False
        if repeat != None:
            if repeat > 0:
                if self._link_type[1] != repeat:
                    return False
            else:
                if self._link_type[1] == -repeat:
                    return False
        if inst != None:
            if inst > 0:
                if self._link_type[2] != inst:
                    return False
            else:
                if self._link_type[2] == -inst:
                    return False
        return True

    def reverse(self):
        assert self.parent()
        l = self.clone()
        l.set_child(self.parent())
        l.set_parent(self.child())
        self.child().del_parent_link(self)
        self.child().add_link(l)
        self.parent().del_link(self)
        self.parent().add_parent_link(l)
        return l

    def copy_data_members_to(self,l):
        l.set_parent(self.parent())
        l.set_parent_pos(self.parent_pos())
        l.set_parent_type(self.parent_type())
        l.set_child(self.child())
        l.set_child_pos(self.child_pos())
        l.set_child_type(self.child_type())
        l.set_id(self.id())
        return l

    def clone(self):
        return self.copy_data_members_to(self.__class__())

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

    def instantiated(self):
        return self._link_type[2] == Linkage.INSTANTIATED

    def is_subst_link(self):
        return self._link_type[0] == Linkage.MONO_TO_SUBST

    def is_subst_out_link(self):
        return self._link_type[0] == Linkage.SUBST_TO_MONO
        
    def is_non_repeat(self):
        return self._link_type[1] == Linkage.NON_REPEAT

    def is_repeat_bridge(self):
        return self._link_type[1] == Linkage.REPEAT_BRIDGE

    def is_repeat_exit(self):
        return self._link_type[1] == Linkage.REPEAT_EXIT

##     def compatible(self,a):
##      if self._parent_type and a._parent_type and self._parent_type != a._parent_type:
##          return False
##      if self._child_type and a._child_type and self._child_type != a._child_type:
##          return False
##      if self._child_pos and a._child_pos and self._child_pos != a._child_pos:
##          return False
##      ppself = self.parentpos2set()
##      ppa    = a.parentpos2set()
##      if len(ppself) > 0 and len(ppa) > 0 and ppself != ppa:
##          return False
##      return True

##     def compatiblewith(self,a):
##      if a._parent_type and self._parent_type != a._parent_type:
##          return False
##      if a._child_type and self._child_type != a._child_type:
##          return False
##      if a._child_pos and self._child_pos != a._child_pos:
##          return False
##      ppself = self.parentpos2set()
##      ppa    = a.parentpos2set()
##      if len(ppa) > 0 and ((len(ppself) == 0) or not (ppself <= ppa)):
##          return False
##      return True

    def equals(self,a):
        # print "---"
        # print repr(self._parent_type), repr(self._child_type), repr(self._child_pos), repr(self._parent_pos)
        # print repr(a._parent_type), repr(a._child_type), repr(a._child_pos), repr(a._parent_pos)
        # print "---"
        # if self._undetermined != a._undetermined:
        #     return False
        if self._link_type != a._link_type:
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
        if not self.instantiated():
            return False
        if self._parent_type == None or len(self._parent_type) != 1:
            return False
        if self._child_type == None or len(self._child_type) != 1:
            return False
        if self._parent_pos == None or len(self._parent_pos) != 1:
            return False
        if self._child_pos == None or len(self._child_pos) != 1:
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
        return (constantString(Linkage,self._link_type[0]),
                constantString(Linkage,self._link_type[1]),
                constantString(Linkage,self._link_type[2]),
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
    _link_type = (Linkage.MONO_TO_SUBST,Linkage.NON_REPEAT,Linkage.INSTANTIATED)

class UninstantiatedLinkage(Linkage):
    _link_type = (Linkage.MONO_TO_MONO,Linkage.NON_REPEAT,Linkage.UNINSTANTIATED)

class SubstOutLinkage(Linkage):
    _link_type = (Linkage.SUBST_TO_MONO,Linkage.NON_REPEAT,Linkage.INSTANTIATED)

class RepeatBridgeLinkage(Linkage):
    _link_type = (Linkage.MONO_TO_MONO,Linkage.REPEAT_BRIDGE,Linkage.INSTANTIATED)

    def __init__(self,**kw):
        super(RepeatBridgeLinkage,self).__init__(self,**kw)
        self._repeat_times_min = 1
        self._repeat_times_max = 1

    def copy_data_members_to(self,l):
        super(RepeatBridgeLink,self).copy_data_members_to(l)
        l.set_repeat_times_min(self.repeat_times_min())
        l.set_repeat_times_max(self.repeat_times_max())
        return l

    # TODO perhaps provide methods for validation the repating times?
    def set_repeat_times_min(self, m):
        assert m is None or isinstance(m, int)
        self._repeat_times_min = m

    def set_repeat_times_max(self, m):
        assert m is None or isinstance(m, int)
        self._repeat_times_max = m

    def repeat_times_min(self):
        return self._repeat_times_min

    def repeat_times_max(self):
        return self._repeat_times_max

class RepeatExitLinkage(Linkage):
    _link_type = (Linkage.MONO_TO_MONO,Linkage.REPEAT_EXIT,Linkage.INSTANTIATED)

class RepeatBridgeSubstOutLinkage(RepeatBridgeLinkage):
    _link_type = (Linkage.SUBST_TO_MONO,Linkage.REPEAT_BRIDGE,Linkage.INSTANTIATED)

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

    print(bdGal)

    bdGlc = bdGal.clone()
    bdGlc.set_id("bdGlc")
    bdGlc.set_stem(Stem.glc)

    print(bdGlc)
    print(bdGal)
    
    #Make GlcNAc...
    GlcNAc = bdGlc.clone()
    GlcNAc.set_id("GlcNAc")
    GlcNAc.add_substituent(Substituent.nAcetyl,
                           parent_pos=2,
                           parent_type=Linkage.oxygenLost,
                           child_type=Linkage.nitrogenAdded)

    print(GlcNAc)
    print(bdGal)
    print(bdGlc)

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
    print(Neu5Ac)

    m = Neu5Ac.clone()
    m.set_id("Clone of Neu5Ac")

    print(m)

    bdGal.add_child(Neu5Ac,
                    parent_type=Linkage.oxygenPreserved,
                    parent_pos=2,
                    child_pos=2,
                    child_type=Linkage.oxygenLost)

    print(bdGal)

