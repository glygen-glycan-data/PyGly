
from . SymbolsTable import SymbolsTable
from . ConstantsTable import ConstantsTable
from . Monosaccharide import *
import re
import sys
import json

from xml.etree.ElementTree import Element, SubElement

class GlycoCTMonoFormat:
    def __init__(self):
        consts = ConstantsTable()
        self.toSym = dict()
        self.fromSym = dict()
        for sym,kv in consts.items():
            if 'GlycoCTSymbol' in kv:
                type,const = sym
                self.toSym[(type,kv['GlycoCTSymbol'])] = const
                self.fromSym[sym] = kv['GlycoCTSymbol']

    def toStr(self,m):
        s = ""
        if m.id() != None:
            s += str(m.id())
        if isinstance(m,Monosaccharide):
            s += "b:"
            s += self.fromSym[('Anomer',m.anomer())]
            if m.stem() != None:
              if m.config() == None:
                cfg = [Config.missing]*len(m.stem())
              else:
                cfg = list(m.config())
              stm = list(m.stem())
              assert len(cfg) == len(stm)
              for cf,st in zip(cfg,stm):
                  s += "-%s%s"%(self.fromSym[('Config',cf)],self.fromSym[('Stem',st)])
            s += "-%s"%self.fromSym[('SuperClass',m.superclass())]
            rs = m.ring_start()
            if rs == None:
                rs = 'x'
            re = m.ring_end()
            if re == None:
                re = 'x'
            s += "-%s:%s"%(rs,re)
            for pi,mi in m.mods():
                s += '|%s:%s'%(",".join(map(str,pi)),self.fromSym[('Mod',mi)])
        elif isinstance(m,Substituent):
            s += "s:"
            s += self.fromSym[('Substituent',m.name())]
        return s

    def toXML(self, m):

        if isinstance(m,Monosaccharide):
            ele = Element('monosaccharide',id=str(m.id()),type="b")
            ele.set('anomer',self.fromSym[('Anomer',m.anomer())])
            ele.set('superclass',self.fromSym[('SuperClass',m.superclass())].lower())
            ele.set('ringStart',str(m.ring_start()) if m.ring_start() != None else "x")
            ele.set('ringEnd',str(m.ring_end()) if m.ring_end() != None else "x")
            ele.set('name',self.toStr(m).split(':',1)[1])
            if m.stem() != None:
                if m.config() == None:
                    cfg = [Config.missing]*len(m.stem())
                else:
                    cfg = list(m.config())
                stm = list(m.stem())
                assert len(cfg) == len(stm)
                for i,(cf,st) in enumerate(zip(cfg,stm)):
                    btele = Element('basetype',id=str(i+1))
                    btele.set('type',self.fromSym[('Config',cf)] + self.fromSym[('Stem',st)])
                    ele.append(btele)
            for pi,mi in m.mods():
                modele = Element('modification',type=self.fromSym[('Mod',mi)])
                modele.set('pos_one',str(pi[0]))
                if len(pi) == 2:
                    modele.set('pos_two',str(pi[1]))
                elif len(pi) > 2:
                    raise RuntimeError("Too many mod positions!")
                ele.append(modele)

        elif isinstance(m,Substituent):
            ele = Element('substituent',id=str(m.id()),type="s",name=self.fromSym[('Substituent',m.name())])

        return ele

    def linkToXML(self, l):
        ele = Element('connection',id=str(l.id()),
                                   parent=str(l.parent().id()),
                                   child=str(l.child().id()))
        linkele = Element('linkage',id=str(l.id()))
        if l.parent_type():
            pt = "|".join(map(lambda t: self.fromSym[('Linkage',t)], sorted(l.parent_type())))
        else:
            pt = self.fromSym[('Linkage',Linkage.missing)]
        linkele.set('parentType',pt)
        if l.child_type():
            ct = "|".join(map(lambda t: self.fromSym[('Linkage',t)], sorted(l.child_type())))
        else:
            ct = self.fromSym[('Linkage',Linkage.missing)]
        linkele.set('childType',ct)
        ele.append(linkele)
        SubElement(linkele,'parent',pos=str(l.posstr(l.parent_pos())))
        SubElement(linkele,'child',pos=str(l.posstr(l.child_pos())))
        return ele

    def linkToStr(self, l, noids=False, replace_parent_id=None, replace_child_id=None):
        # 1:1d(2+1)2n
        s = ""
        if l.id() != None and not noids:
            s += str(l.id())+":"
        if l.parent().id() != None and not noids:
            if replace_parent_id != None:
                s += str(replace_parent_id)
            else:
                s += str(l.parent().id())
        if l.parent_type():
            s += "|".join(map(lambda t: self.fromSym[('Linkage',t)], sorted(l.parent_type())))
        else:
            s += self.fromSym[('Linkage',Linkage.missing)]
        s += '(%s'%(l.posstr(l.parent_pos()),)
        s += '+%s)'%(l.posstr(l.child_pos()),)
        if l.child().id() != None and not noids:
            if replace_child_id != None:
                s += str(replace_child_id)
            else:
                s += str(l.child().id())
        if l.child_type():
            s += "|".join(map(lambda t: self.fromSym[('Linkage',t)], sorted(l.child_type())))
        else:
            s += self.fromSym[('Linkage',Linkage.missing)]
        return s

    fromStrRE = re.compile(r'^(\d+)([bs]):(.*)$')
    def fromStr(self,mstr):
        m = self.fromStrRE.search(mstr.strip())
        if not m:
            raise RuntimeError("Bad GlycoCT monosaccharide line:"+mstr)
        id = int(m.group(1))
        type = m.group(2)
        desc = m.group(3)

        if type == 's':

            s = Substituent(self.toSym[('Substituent',desc)])
            s.set_id(id)
            return s

        m = Monosaccharide()
        m.set_id(id)
        m.set_external_descriptor_id(id)
        MODS = desc.split('|')
        MONO = MODS.pop(0).split('-')
        
        #set ring
        ring =  MONO.pop()
        ringStart,ringEnd = ring.split(':')
        try:
            m.set_ring_start(int(ringStart))
        except ValueError:
            if ringStart != "x":
                raise RuntimeError('Bad GlycoCT monosaccharide line - bad ring start value: '+mstr)
        try:
            m.set_ring_end(int(ringEnd))
        except ValueError:
            if ringEnd != "x":
                raise RuntimeError('Bad GlycoCT monosaccharide line - bad ring end value: '+mstr)

        #set class
        superclass = MONO.pop()
        try:
            m.set_superclass(self.toSym[('SuperClass',superclass)])
        except KeyError:
            raise RuntimeError('Bad GlycoCT monosaccharide line - unsupported superclass: '+mstr)

        #set anomer
        anomer = MONO.pop(0)
        m.set_anomer(self.toSym[('Anomer',anomer)])

        #set stem
        configs = []
        stems = []
        for st in MONO:
            configs.append(self.toSym[('Config',st[0])])
            stems.append(self.toSym[('Stem',st[1:])])
        
        if len(stems) > 0:
            if zip(configs,stems) != [(None,None)]:
                m.set_config(*configs)
                m.set_stem(*stems)

        #set mods
        for mod in MODS:
            num,modi = mod.split(':')
            modi=self.toSym[('Mod',modi)]
            m.add_mod(num,modi)
        return m

    linkFromStrRE = re.compile(r'^(REP)?(\d+):(\d+)([dnohx])[(](.+)\+(.+)[)](\d+)([dnohx])(=(-?\d+)-(-?\d+))?$')
    def linkParseBase(self, s):
        m = self.linkFromStrRE.search(s)
        if not m:
            raise RuntimeError("Bad GlycoCT link line:" + s)

        repeatbridge = m.group(1) == "REP"

        id = int(m.group(2))
        parentid = int(m.group(3))
        parenttype = self.toSym[('Linkage', m.group(4))]
        parentpos2 = None
        try:
            parentpos = map(int, m.group(5).split('|'))
            if -1 in parentpos:
                parentpos = None
        except ValueError:
            parentpos = None
        try:
            childpos = int(m.group(6))
            if childpos < 1:
                childpos = None
        except ValueError:
            childpos = None
        childid = int(m.group(7))
        childtype = self.toSym[('Linkage', m.group(8))]

        repeattimemin, repeattimemax = None, None
        if repeatbridge:
            repeattimemin = int(m.group(10))
            repeattimemax = int(m.group(11))

            if repeattimemin == -1:
                repeattimemin = None
            if repeattimemax == -1:
                repeattimemax = None

        return id, parentid, parentpos, parenttype, childid, childpos, childtype, repeatbridge, repeattimemin, repeattimemax

    def linkFromPara(self, parent, child, parenttype, parentpos, childtype, childpos, repeat_exit=False, repeat_bridge=False, undet=False):
        assert not (repeat_exit and repeat_bridge)
        link_characteristics = dict(parent_type=parenttype, parent_pos=parentpos, 
                                    child_type=childtype, child_pos=childpos)
        if child.is_monosaccharide():
            if repeat_exit:
                linkage = RepeatExitLinkage(**link_characteristics)
                l = parent.add_child_with_special_linkage(child, linkage)
            elif repeat_bridge:
                if parent.is_substituent():
                    linkage = RepeatBridgeSubstOutLinkage(**link_characteristics)
                else:
                    linkage = RepeatBridgeLinkage(**link_characteristics)
                l = parent.add_child_with_special_linkage(child, linkage)
            elif undet:
                linkage = UninstantiatedLinkage(**link_characteristics)
                l = parent.add_child_with_special_linkage(child, linkage)
            elif parent.is_substituent():
                linkage = SubstOutLinkage(**link_characteristics)
                l = parent.add_child_with_special_linkage(child, linkage)
            else:
                l = parent.add_child(child, **link_characteristics)
        else:
            if child.name() == Substituent.amino and parenttype == Linkage.oxygenPreserved:
                child._sub = Substituent.amino_oxygen_preserved
            elif child.name() == Substituent.methyl and parenttype == Linkage.oxygenLost:
                child._sub = Substituent.methyl_oxygen_lost
            elif child.name() == Substituent.phosphate:
                if parenttype == Linkage.oxygenLost:
                    child._sub = Substituent.phosphate_oxygen_lost
                elif child.has_parent_links():
                    # already an edge to this substituent...
                    child._sub = Substituent.phosphate_bridged
            elif child.name() == Substituent.sulfate and parenttype == Linkage.oxygenLost:
                child._sub = Substituent.sulfate_oxygen_lost
            elif child.name() == Substituent.acetyl and parenttype == Linkage.oxygenLost:
                child._sub = Substituent.acetyl_oxygen_lost
            l = parent.add_substituent(child,**link_characteristics)
        return l

    def linkFromStr(self, s, res, undet=False):

        id, parentid, parentpos, parenttype, childid, childpos, childtype, r1, r2, r3 = self.linkParseBase(s)

        if parentid not in res:
            raise RuntimeError("Bad GlycoCT link line, parent missing:"+s)
        if childid not in res:
            raise RuntimeError("Bad GlycoCT link line, child missing:"+s)
        if parentid >= childid:
            raise RuntimeError("Bad GlycoCT link line, backwards link:"+s)
        parent = res[parentid]
        #if isinstance(parent, Substituent):
        #    raise RuntimeError("Bad GlycoCT link line, substituent as parent:"+s)
        child = res[childid]

        l = self.linkFromPara(parent, child, parenttype, parentpos, childtype, childpos, undet=undet)

        return [l]

class MonosaccharideDB(GlycoCTMonoFormat):
   def toStr(self,m):
        s = [ GlycoCTMonoFormat.toStr(self,m).split(':',1)[1] ]
        for l in m.substituent_links():
            sub = l.child()
            ls0 = self.linkToStr(l,noids=True)
            f,t = ls0.split('+')
            ls = ''.join(f.split('(')[::-1]) + ":" + ''.join(t.split(')'))
            # print(ls)
            # ls = "%s%s:%s%s"%(list(l.parent_pos())[0],self.fromSym[('Linkage',list(l.parent_type())[0])],list(l.child_pos())[0],self.fromSym[('Linkage',list(l.child_type())[0])])
            # print(ls)
            ss = GlycoCTMonoFormat.toStr(self,sub).split(':',1)[1]
            s.append('(%s)%s'%(ls,ss))
        return "||".join(s)

class MonoSymLookup(dict):
    def __init__(self):
        st = SymbolsTable()
        for key,kv in st.items():
            if self.__class__.__name__ in kv:
                self[key] = kv[self.__class__.__name__]
    def key(self,m):
        if isinstance(m,Monosaccharide):
            supcls = tuple(sorted(('SuperClass',s) for s in [m.superclass()]) if m.superclass() != None else ())
            stem = tuple(sorted(('Stem',s) for s in m.stem()) if m.stem() != None else ())
            mods = tuple(sorted(('Mod',m[1]) for m in m.mods()) if m.mods() != None else ())
            subst = tuple(sorted(('Substituent',s.name()) for s in m.substituents()))
        elif isinstance(m,Substituent):
            supcls = ()
            stem = ()
            mods = ()
            subst = (('Substituent',m.name()),)
        return supcls,stem,mods,subst
    def toStr(self,m):
        if m != None:
            k = self.key(m)
        else:
            raise RuntimeError("Monosaccharide is None.")
        try:
            return self[k]
        except KeyError:
            # if self.__class__.__name__ == 'IUPACSym':
            #     print >>sys.stderr, "%s table can't find: %s"%(self.__class__.__name__,k)
            #     print >>sys.stderr, '\n'.join(map(str,sorted(self.keys())))
            # sys.exit(1)
            raise

class IUPACSym(MonoSymLookup):
    pass

class LinCodeSym(MonoSymLookup):
    pass

class LinCodeRank(MonoSymLookup):
    def __init__(self):
        st = SymbolsTable()
        for key,kv in st.items():
            if self.__class__.__name__ in kv:
                self[key] = int(kv[self.__class__.__name__])

class MassSym(MonoSymLookup):
    pass

class GlycamSym(MonoSymLookup):
    pass

class JSONCanonResidue(object):
    iupac = IUPACSym()
    mdb = MonosaccharideDB()
    gctmf = GlycoCTMonoFormat()
    ringToChar = {(1,5): "p", (2,6): "p", (0,0): "ol", (1,4): "f"}
    configToChar = {"d": "D", "l": "L", ("d","d"): "D"}
    def toStr(self,m):
        data = dict()
        if isinstance(m,Monosaccharide):
            data['residueid'] = str(m.id())
            data['residuetype'] = 'monosaccharide'
            data['monodb'] = self.mdb.toStr(m)
            try:
                data['name'] = self.iupac.toStr(m)
            except KeyError:
                pass #data['iupac'] = 'Xxx'
            if m.anomer() != Anomer.missing:
                data['anomer'] = self.gctmf.fromSym[('Anomer',m.anomer())]
            else:
                pass #data['anomer'] = 'x'
            # if m.superclass() != SuperClass.missing:
            #     data['superclass'] = self.gctmf.fromSym[('SuperClass',m.superclass())].lower()
            # else:
            #     pass #data['superclass'] = 'xxx'
            rs = m.ring_start()
            re = m.ring_end()
            ringkey = (rs,re)
            if ringkey in self.ringToChar:
                data['ring'] = self.ringToChar[ringkey]
            cfst = []
            if m.stem() != None:
                if m.config() == None:
                    cfg = [Config.missing]*len(m.stem())
                else:
                    cfg = list(m.config())
                stm = list(m.stem())
                assert len(cfg) == len(stm)
                
                for i,(cf,st) in enumerate(zip(cfg,stm)):
                    cfst.append((self.gctmf.fromSym[('Config',cf)],self.gctmf.fromSym[('Stem',st)]))
                
            if len(cfst) == 1 and cfst[0][0] in self.configToChar:
                data['absolute'] = self.configToChar[cfst[0][0]]
            if len(cfst) == 2 and (cfst[0][0],cfst[1][0]) in self.configToChar and data.get('iupac') in ("NeuAc","NeuGc"):
                data['absolute'] = self.configToChar[(cfst[0][0],cfst[1][0])]

            if m.parent_link_count() == 1:
                pl = m.any_parent_link()
                if pl.parent_pos() != None:
                    data['site'] = str(pl.posstr(pl.parent_pos()))
                data['parentid'] = str(pl.parent().id())
      
        elif isinstance(m,Substituent):

            data['residueid'] = str(m.id())
            data['residuetype'] = 'substituent'
            data['monodb'] = self.mdb.toStr(m)
            data['name'] = self.gctmf.toStr(m).split(':',1)[1]

        else:
            raise NotImplemented()

        return json.dumps(data,indent=2,sort_keys=True);

