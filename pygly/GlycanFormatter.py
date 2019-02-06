
from MonoFormatter import GlycoCTMonoFormat, LinCodeSym, LinCodeRank, IUPACSym, GlycamSym
from Monosaccharide import Monosaccharide, Linkage, Anomer, Substituent, Mod
from Glycan import Glycan
from MonoFactory import MonoFactory

import re, sys
import copy
from collections import defaultdict

class GlycanFormatter:
    def writeToFile(self,thefile,glycan):
        doclose = False
        if isinstance(thefile,basestring):
            thefile = open(thefile,'w')
            doclose = True
        thefile.write(self.toStr(glycan))
        if doclose:
            thefile.close()
    def readFromFile(self,thefile):
        pass

class GlycanParseError(Exception):
    def __str__(self):
        return self.message

class GlycoCTParseError(GlycanParseError):
    pass

class GlycoCTRepeatedSectionError(GlycoCTParseError):
    def __init__(self,section,lineno):
	self.section = section
	self.message = "GlycoCT parser, line %d: Repeated section %s"%(lineno, section)

class GlycoCTLINBeforeRESError(GlycoCTParseError):
    def __init__(self,lineno):
	self.message = "GlycoCT parser, line %d: LIN section before RES section"%(lineno,)

class GlycoCTUNDBeforeRESLINError(GlycoCTParseError):
    def __init__(self,lineno):
	self.message = "GlycoCT parser, line %d: UND section before RES or LIN section"%(lineno,)

class GlycoCTUnsupportedSectionError(GlycoCTParseError):
    def __init__(self,section,lineno):
	self.section = section
	self.message = "GlycoCT parser, line %d: Unsupported section %s"%(lineno, section)

class GycoCTBadRESLineError(GlycoCTParseError):
    def __init__(self,message,lineno,line):
	self.message = "GlycoCT parser, line %d: %s"%(lineno, message)
	self.line = line

class GlycoCTBadLINLineError(GlycoCTParseError):
    def __init__(self,message,lineno,line):
	self.message = "GlycoCT parser, line %d: %s"%(lineno, message)
	self.line = line

class GlycoCTParentLinkError(GlycoCTParseError):
    def __init__(self,message,lineno,line):
	self.message = "GlycoCT parser, line %d: %s"%(lineno, message)
	self.line = line

class GlycoCTUnexpectedLineError(GlycoCTParseError):
    def __init__(self,lineno,line):
	self.message = "GlycoCT parser, line %d: Unexpected line:\n ==> %s"%(lineno,line)
	self.line = line

class GlycoCTUndeterminedLinkageError(GlycoCTParseError):
    def __init__(self):
	self.message = "GlycoCT parser, undetermined linkage instantiation error"

class GlycoCTUnconnectedCountError(GlycoCTParseError):
    def __init__(self):
	self.message = "GlycoCT parser, unconnected node count error"

class GlycoCTFormat(GlycanFormatter):
    def __init__(self):
        self.monofmt = GlycoCTMonoFormat()
    def toStr(self,g):

        g.set_ids()
        
        roots = []
        if g.root():
            roots.append(g.root())
        for ur in g.undetermined_roots():
	    if not ur.connected():
                roots.append(ur)

        s = "RES\n"
        for r in roots:
            if len(r.parent_links()) == 0 or r == g.root():
                for m in g.subtree_nodes(r,subst=True):
                    s += self.monofmt.toStr(m)+"\n"
        
        first = True
        linkid = 1
        for r in roots:
            if len(r.parent_links()) == 0 or r == g.root():
                for l in sorted(g.subtree_links(r,subst=True),key=lambda l: l.child().id()):
                    if first:
                        s += "LIN\n"
                        first = False
                    l.set_id(linkid)
                    linkid += 1
                    s += self.monofmt.linkToStr(l)+"\n"

        first = True
        undind = 0
        for r in roots:
            # print r, g.root(), r.parent_links()
            if len(r.parent_links()) == 0:
                continue
            if r == g.root():
                continue
            if first:
                s += "UND\n"
		first = False
            undind += 1
            s += "UND%s:100.0:100.0\n"%(undind,)
            parentids = set()
            subtreelinkage = set()
            for pl in r.parent_links():
                parentids.add(pl.parent().id())
                subtreelinkage.add(self.monofmt.linkToStr(pl,noids=True))
            s += "ParentIDs:%s\n"%("|".join(map(str,sorted(parentids))),)
            assert len(subtreelinkage) == 1
            s += "SubtreeLinkageID1:%s\n"%(subtreelinkage.pop(),)
            s += "RES\n"
            for m in g.subtree_nodes(r,subst=True):
                s += self.monofmt.toStr(m)+"\n"
            first1 = True
            for l in sorted(g.subtree_links(r,subst=True),key=lambda l: l.child().id()):
                if first1:
                    s += "LIN\n"
                    first1 = False
                l.set_id(linkid)
                linkid += 1
                s += self.monofmt.linkToStr(l)+"\n"
        return s

    def toGlycan(self,s):
        res = {}
	und = defaultdict(dict)
	undind = None
	state = None
        undets = set()
	seen = set()
        for lineno,l in enumerate(s.splitlines()):
            l = l.strip()
	    if not l:
		continue
            if l == "RES":
		if l in seen and state != "UND":
		    raise GlycoCTRepeatedSectionError(lineno=lineno+1,section=l)
		if state == "UND":
                    state = "UNDRES"
		else:
		    state = "RES"
		seen.add(state)
                continue
            if l == "LIN":
		if l in seen and state != "UNDRES":
		    raise GlycoCTRepeatedSectionError(lineno=lineno+1,section=l)
		if "RES" not in seen:
		    raise GlycoCTLINBeforeRESError()
		if state == "UND":
		    state = "UNDLIN"
		else:
		    state = "LIN"
		seen.add(state)
                continue
	    if l == "UND":
		if l in seen:
		    raise GlycoCTRepeatedSectionError(lineno=lineno+1,section=l)
		if "RES" not in seen or "LIN" not in seen:
		    raise GlycoCTUNDBeforeRESLINError()
		state = "UND"
		seen.add(state)
		continue
	    if re.search(r'^UND\d+:',l):
		m = re.search(r"^UND(\d+):",l)
		state = "UND"
		if m:
		    undind = int(m.group(1))	
		    und[undind]['frac'] = map(float,l.split(':')[1:])
		    continue
	    if re.search(r'^[A-Z][A-Z][A-Z]$',l):
		raise GlycoCTUnsupportedSectionError(lineno=lineno+1,section=l)
            if state in ("RES","UNDRES"):
		try:
                    m = self.monofmt.fromStr(l)
		except RuntimeError, e:
		    raise GycoCTBadRESLineError(message=e.message,lineno=lineno+1,line=l)
                res[m.id()] = m
		if state == "UNDRES" and undind != None and 'root' not in und[undind]:
		    und[undind]['root'] = m.id()
                    undets.add(m)
                continue
            if state in ("LIN","UNDLIN"):
		try:
                    links = self.monofmt.linkFromStr(l,res)
		    if len(links) > 1:
			for l in links:
			    l.set_undetermined(True)
			    l.set_instantiated(False)
			    l.child().set_connected(False)
		except (RuntimeError,AttributeError), e:
		    raise GlycoCTBadLINLineError(message=e.message,lineno=lineno+1,line=l)	
                continue
	    if state == "UND":
		m = re.search(r"^UND(\d+):",l)
		if m:
		    undind = int(m.group(1))	
		    und[undind]['frac'] = map(float,l.split(':')[1:])
		    continue
		m = re.search(r"^ParentIDs:(\d+(\|\d+)*)$",l)
		if m:
		    und[undind]['parentids'] = map(int,m.group(1).split('|'))
		    continue
		m = re.search(r"^SubtreeLinkageID(\d+):(.*)$",l)
		if m:
		    und[undind]['stlink'] = (m.group(1),m.group(2))
		    continue
	    raise GlycoCTUnexpectedLineError(lineno=lineno+1,line=l)
	for d in und.values():
	    linkagestr = d['stlink'][1]
	    rootid = d['root']
	    for pid in d['parentids']:
		try:
		  links = self.monofmt.linkFromStr("0:%s%s%s%s"%(pid,linkagestr[:-1],rootid,linkagestr[-1]),res)
		  for l in links:
		    l.set_undetermined(True)
		    l.set_instantiated(False)
		    l.child().set_connected(False)
		except (RuntimeError,AttributeError), e:
		    raise GlycoCTParentLinkError(message=e.message,lineno=lineno+1,line=l)	
        unconnected = set()
        monocnt = 0
        for id,r in res.items():
            if not isinstance(r,Monosaccharide):
                continue
            monocnt += 1
            if len(r.parent_links()) > 0:
                continue
            r.set_connected(False)
            unconnected.add(r)
	if len(unconnected) not in (1,monocnt):
	    raise GlycoCTUnconnectedCountError()
        # single monosacharides are considered a structure, not a composition...
        if len(unconnected) == 1:
            g = Glycan(unconnected.pop())
            g.set_undetermined(undets)
        else:
            assert len(undets) == 0
            g = Glycan()
            g.set_undetermined(unconnected)
        return g

class LinearCodeParseError(GlycanParseError):
    pass

class LinearCodeBadFormat(LinearCodeParseError):
    def __init__(self,code,pos):
	self.code = code
	self.pos = pos
	self.message = "Linear Code parser, position %d: Bad format in linear code: %s"%(pos,code)

class LinearCodeBadSym(LinearCodeParseError):
    def __init__(self,code,pos,sym):
	self.code = code
	self.pos = pos
	self.sym = sym
	self.message = "Linear Code parser, position %d: Bad symbol %s linear code: %s"%(pos,sym,code)

class LinearCodeBadAnomer(LinearCodeParseError):
    def __init__(self,code,pos,anomer):
	self.code = code
	self.pos = pos
	self.anomer = anomer
	self.message = "Linear Code parser, position %d: Bad anomer %s in linear code: %s"%(pos,anomer,code)

class LinearCodeBadPosition(LinearCodeParseError):
    def __init__(self,code,pos,badpos):
	self.code = code
	self.pos = pos
	self.badpos = badpos
	self.message = "Linear Code parser, position %d: Bad position %s in linear code: %s"%(pos,badpos,code)

class LinearCodeFormat(GlycanFormatter):
    def __init__(self,connection=None):
        self.monofmt = LinCodeSym()
        self.lcrank = LinCodeRank()
        self.mf = MonoFactory()
    def toStr(self,g):
        return self.linearstr(None,g.root())
    def linchkey(self,l):
        pos = l.parent_pos()
        if not pos:
            pos = -100
        rank = self.lcrank.toStr(l.child())
        return (rank,-pos)
    def linearstr(self,pl,m):
        s = self.monofmt.toStr(m)
        if not pl:
            conn = ""
        else:
            pos=pl.parent_pos()
            if not pos:
                pos = '?'
            if m.anomer() == Anomer.alpha:
                an = 'a'
            elif m.anomer() == Anomer.beta:
                an = 'b'
            elif m.anomer() == Anomer.uncyclized:
                an = 'o'
            else:
                an = '?'
            conn = (an + str(pos))
        kidlinks = sorted(m.links(),key=self.linchkey)
        n = len(kidlinks)
        if n == 0:
            return s+conn
        elif n == 1:
            return self.linearstr(kidlinks[0],kidlinks[0].child())+s+conn
        else:
            t = ""
            for cl in kidlinks[:-1]:
                t = "(" + self.linearstr(cl,cl.child()) + ")" + t
            return self.linearstr(kidlinks[-1],kidlinks[-1].child())+t+s+conn

    def toGlycan(self, s):
	orig = s
 	m = re.search(r'([A-Z]+)$',s)
	if not m:
	    raise LinearCodeBadFormat(code=orig,pos=len(s))
	sym = m.group(1)
        lcsym = "LinearCode:"+sym
        if lcsym in self.mf:
            root = self.mf.new(lcsym)
        else:
            raise LinearCodeBadSym(code=orig,pos=len(s),sym=sym)
	s = s[:-len(sym)]

        branchpoint = []
        parent = root
        while s != "":
            # either a number or a bracket.
            if s[-1] not in '()':
	        m = re.search(r'([A-Z]+)(.)(.)$',s)
		if not m:
		    raise LinearCodeBadFormat(code=orig,pos=len(s))
		sym = m.group(1)
		lcsym = "LinearCode:"+sym
                anomer = m.group(2)
		pos = m.group(3)
                if lcsym not in self.mf:
                    raise LinearCodeBadSym(code=orig,pos=len(s),sym=sym)
                if anomer not in 'ab?':
                    raise LinearCodeBadAnomer(code=orig,pos=len(s),anomer=anomer)
                if pos not in '?123456789':
                    raise LinearCodeBadPosition(code=orig,pos=len(s),badpos=pos)
                if anomer == 'a':
                    anomer = Anomer.alpha
                elif anomer == 'b':
                    anomer = Anomer.beta
                else:
                    anomer = Anomer.missing
                if pos == '?':
                    pos = None
                else:
                    pos = int(pos)
		m = self.mf.new(lcsym)
                m.set_anomer(anomer)
                parent.add_child(m,parent_pos=pos)
		s = s[:-(2+len(sym))]
                parent = m
            elif s[-1] == ')':
                branchpoint.append(parent)
                s=s[:-1]
            elif s[-1] == '(':
                parent = branchpoint.pop()
                s=s[:-1]
            else:
                raise RuntimeError("Bad linear code format!")
        return Glycan(root)



class IUPACLinearParseError(GlycanParseError):
    pass

class IUPACLinearBadFormat(IUPACLinearParseError):
    def __init__(self,code,pos):
	self.code = code
	self.pos = pos
	self.message = "IUPAC Linear parser, position %d: Bad format in linear code: %s^%s."%(pos,code[:pos],code[pos:])

class IUPACLinearBadSym(IUPACLinearParseError):
    def __init__(self,code,pos,sym):
	self.code = code
	self.pos = pos
	self.sym = sym
	self.message = "IUPAC Linear parser, position %d: Bad symbol %s linear code: %s"%(pos,sym,code)

class IUPACLinearBadAnomer(IUPACLinearParseError):
    def __init__(self,code,pos,anomer):
	self.code = code
	self.pos = pos
	self.anomer = anomer
	self.message = "IUPAC Linear parser, position %d: Bad anomer %s in linear code: %s"%(pos,anomer,code)

class IUPACLinearBadPosition(IUPACLinearParseError):
    def __init__(self,code,pos,badpos):
	self.code = code
	self.pos = pos
	self.badpos = badpos
	self.message = "IUPAC Linear parser, position %d: Bad position %s in linear code: %s"%(pos,badpos,code)

class IUPACLinearFormat(GlycanFormatter):
    def __init__(self,connection=None):
        self.monofmt = IUPACSym()
        self.mf = MonoFactory()
    def toStr(self,g):
        return self.linearstr(None,g.root())
    def linkkey(self,l):
        pos = l.parent_pos()
        if not pos:
            pos = -100
        return -pos
    def linearstr(self,pl,m):
        aglycon = None
        if not pl:
            try:
                s = self.monofmt.toStr(m)
            except KeyError:
                if m.noring():
                    if m.count_mod() == 1 and \
                       m.count_mod(Mod.aldi) == 1:
                        aglycon = 'aldi?'
                        m1 = m.clone()
                        m1.clear_mods()
                        s = self.monofmt.toStr(m1)
                    else:
                        raise
                else:
                    raise
        else:
            s = self.monofmt.toStr(m)
        if not pl:
	    if m.anomer() == Anomer.alpha:
		conn = 'a'
	    elif m.anomer() == Anomer.beta:
		conn = 'b'
	    else:
		conn = ""
        else:
            pos=pl.parent_pos()
            if not pos:
                pos = '?'
            cpos = pl.child_pos()
            if not pos:
                cpos = '?'
            if m.anomer() == Anomer.alpha:
                an = 'a'
            elif m.anomer() == Anomer.beta:
                an = 'b'
            else:
                an = '?'
            conn = (an + str(cpos) + "-" + str(pos))
        kidlinks = sorted(m.links(),key=self.linkkey)
        n = len(kidlinks)
        if n == 0:
            return s+conn
        elif n == 1:
            return self.linearstr(kidlinks[0],kidlinks[0].child())+s+conn
        else:
            t = ""
            for cl in kidlinks[:-1]:
                t = "(" + self.linearstr(cl,cl.child()) + ")" + t
            return self.linearstr(kidlinks[-1],kidlinks[-1].child())+t+s+conn

    def toGlycan(self, s):
	orig = s
	m = re.search(r'(?P<subst>(\d[SP])?\(\d[SP]\)+)?(?P<sym>(Glc|Gal|Man|Fuc|Xyl)[a-zA-Z5926,]+)$',s)
	if not m:
	    raise IUPACLinearBadFormat(code=orig,pos=len(s))
	sym = m.group('sym')
        if sym in self.mf:
	    anomer = None
            root = self.mf.new(sym)
	elif sym[-1] in 'ab' and sym[:-1] in self.mf:
	    anomer = sym[-1]
	    root = self.mf.new(sym[:-1])
        else:
            raise IUPACLinearBadSym(code=orig,pos=len(s),sym=sym)
	if anomer:
	     if anomer == 'a':
                  anomer = Anomer.alpha
             elif anomer == 'b':
                  anomer = Anomer.beta
	     root.set_anomer(anomer)
        if m.group('subst'):
            for subst in re.findall(r'\d[SP]',m.group('subst')):
                if subst[1] == 'S':
                    root.add_substituent(Substituent.sulfate,
                                         parent_pos=int(subst[0]),
					 child_pos=1,
					 parent_type=Linkage.oxygenPreserved,
					 child_type=Linkage.nitrogenAdded)
                elif subst[1] == 'P':
                    root.add_substituent(Substituent.phosphate,
                                         parent_pos=int(subst[0]),
					 child_pos=1,
					 parent_type=Linkage.oxygenPreserved,
					 child_type=Linkage.nitrogenAdded)
	s = s[:-len(m.group(0))]

        branchpoint = []
        parent = root
        while s != "":
            # either a number or a bracket.
            if s[-1] not in '()':
	        m = re.search(r'(?P<subst>(\d[SP])?(\(\d[SP]\))+)?(?P<sym>[A-Z][a-zA-Z592,]+)(?P<anomer>.)(?P<cpos>.)-(?P<ppos>.)$',s)
                if not m:
		    raise IUPACLinearBadFormat(code=orig,pos=len(s))
		sym = m.group('sym')
		anomer = m.group('anomer')
		cpos = m.group('cpos')
                ppos = m.group('ppos')
                subst = m.group('subst')
                remove = len(m.group(0))
                if sym not in self.mf:
                    raise IUPACLinearBadSym(code=orig,pos=len(s),sym=sym)
                if anomer not in 'ab?':
                    raise IUPACLinearBadAnomer(code=orig,pos=len(s),anomer=anomer)
                if cpos not in '?123456789':
                    raise IUPACLinearBadPosition(code=orig,pos=len(s),badpos=cpos)
                if ppos not in '?123456789':
                    raise IUPACLinearBadPosition(code=orig,pos=len(s),badpos=ppos)
                if anomer == 'a':
                    anomer = Anomer.alpha
                elif anomer == 'b':
                    anomer = Anomer.beta
                else:
                    anomer = Anomer.missing
                if cpos == '?':
                    cpos = None
                elif cpos:
                    cpos = int(cpos)
                if ppos == '?':
                    ppos = None
                else:
                    ppos = int(ppos)
		m = self.mf.new(sym)
                m.set_anomer(anomer)
                parent.add_child(m,child_pos=cpos,parent_pos=ppos,
				   parent_type=Linkage.oxygenPreserved,
                                   child_type=Linkage.oxygenLost)
                if subst:
                    for si in re.findall(r'\d[SP]',subst):
                        if si[1] == 'S':
                            m.add_substituent(Substituent.sulfate,
                                              parent_pos=int(si[0]),
					      child_pos=1,
					      parent_type=Linkage.oxygenPreserved,
					      child_type=Linkage.nitrogenAdded)
                        elif si[1] == 'P':
                            m.add_substituent(Substituent.phosphate,
                                              parent_pos=int(si[0]),
					      child_pos=1,
					      parent_type=Linkage.oxygenPreserved,
					      child_type=Linkage.nitrogenAdded)
		s = s[:-remove]
                parent = m

            elif s[-1] == ')':
                branchpoint.append(parent)
                s=s[:-1]
            elif s[-1] == '(':
                parent = branchpoint.pop()
                s=s[:-1]
            else:
                raise RuntimeError("Bad IUPAC linear format!")
        return Glycan(root)

class IUPACParserAbstract():
    # Link tuple (parent_pos, child_pos)
    regexResTemplate = {
        "anomer": None,
        "skel": None,
        "link": tuple(),
        "substituent": None,
        "sublink": 0,
        "branchpointstart": 0,
        "branchpointend": 0,
    }

    description = "This is the abstract class for IUPAC parser"
    example = ""  # An example IUPAC string
    alias = ""  # Used as prefix to extract monosaccharide info from ini.

    precompiledpattern = r""  # Need to be overwritten by each IUPAC Parser
    pattern = ""  # Compiled pattern

    def __init__(self):
        self.mf = MonoFactory()
        self.pattern = re.compile(self.precompiledpattern)

    def toGlycan(self, seq):
        # Preparation "Global Variable of this method"
        root = 0
        branchPoint = []
        parentMono = 0

        # Process
        regexRes = self.regexSearch(seq)
        for i, monoinfo in enumerate(regexRes):
            mono = self.monoAssemble(monoinfo)
            add_stack = monoinfo["branchpointstart"]
            remove_stack = monoinfo["branchpointend"]

            if i == 0:
                root = mono
                parentMono = mono

            else:
                # Debug: add parent_pos, child_pos and link_type
                cp, pp = monoinfo["link"]
                if pp:
                    try:
                        pp = map(lambda x: int(x), pp.split("/"))
                    except ValueError:
                        raise IUPACLinearBadPosition(code=seq, pos=len(seq), badpos=pp)
                if cp:
                    try:
                        cp = map(lambda x: int(x), cp.split("/"))
                    except ValueError:
                        raise IUPACLinearBadPosition(code=seq, pos=len(seq), badpos=cp)

                parentMono.add_child(mono, child_pos=cp, parent_pos=pp,
                                     parent_type=Linkage.oxygenPreserved,
                                     child_type=Linkage.oxygenLost)

                # Debug: Possibility for multiple branching starting bracket
                # It doesn't really make any sense to be honest.
                if add_stack != 0:
                    if add_stack > 1:
                        raise IUPACLinearBadPosition("Multiple branching starting bracket", 0, 0)
                    branchPoint.append(parentMono)
                parentMono = mono

                # Well... It can make sense if there are multiple ending bracket. Just pop n element from branchPoint
                if remove_stack != 0:
                    if remove_stack > 1:
                        raise IUPACLinearBadPosition("Multiple branching ending bracket", 0, 0)
                    parentMono = branchPoint.pop(-1)
        glycanObj = Glycan(root)
        return glycanObj

    def patternCompilation(self):
        raise NotImplemented

    def regexSearch(self, seq):
        raise NotImplemented

    def monoAssemble(self, info):
        fullname = self.alias + info["skel"]
        try:
            m = self.mf.new(fullname)
        except:
            raise LinearCodeBadSym("", 0, fullname)
        if m.anomer() != Anomer.uncyclized and info["anomer"]:
            m.set_anomer(eval(info["anomer"]))
        return m


class IUPACParserGlyTouCanCondensed(IUPACParserAbstract):
    description = "Usually used by GlyTouCan.org"
    example = "Man(a1-3/6)[Man(a1/3/5-?)Gal(b1-3)]GalNAc(?1-"

    alias = "GTCC:"
    precompiledpattern = r"(?P<bpe>(\[)*)(?P<mono>(Glc|Gal|Man|Fuc|Xyl)([a-zA-Z5926,]+)?)(?P<link>(\([ab?]?[1-9?](/[1-9])*-([1-9?](/[1-9])*\))?))(?P<bps>(\])*)"
    pattern = ""

    def patternCompilation(self):
        # Deprecated
        pm = r'(?P<mono>(Glc|Gal|Man|Fuc|Xyl)([a-zA-Z5926,]+)?)'
        pl = r'(?P<link>(\([ab?]?[1-9?](/[1-9])*-([1-9?](/[1-9])*\))?))'
        pc = r"(?P<bpe>(\[)*)" + pm + pl + r"(?P<bps>(\])*)"
        compiledpattern = re.compile(pc)
        self.pattern = compiledpattern

    def regexSearch(self, seq):
        searchres = [m.groupdict() for m in self.pattern.finditer(seq)]
        searchres = list(reversed(searchres))

        res = []
        for s in searchres:
            res0 = copy.deepcopy(self.regexResTemplate)

            anomersl = s["link"][1]
            if anomersl == "a":
                anomer = "Anomer.alpha"
            elif anomersl == "b":
                anomer = "Anomer.beta"
            elif anomersl == "?":
                anomer = None
            else:
                # Unsupported anomer
                raise ValueError

            link = s["link"].lstrip("(").rstrip(")")
            link = link[1:].split("-")
            link = map(lambda x: None if x == "" else x, link)
            link = map(lambda x: None if x == "?" else x, link)

            res0["anomer"] = anomer
            res0["skel"] = s["mono"]
            res0["link"] = tuple(link)
            res0["substituent"] = None
            res0["sublink"] = None
            res0["branchpointstart"] = len(s["bps"])
            res0["branchpointend"] = len(s["bpe"])
        return searchres


class IUPACParserGlyTouCanExtended(IUPACParserAbstract):
    description = "Usually used by GlyTouCan.org"
    example = "alpha-D-Galp-(1->3)-D-Gal(?->"

    alias = "GTCE:"
    precompiledpattern = r"(?P<bpe>(\[)*)(?P<anomer>(alpha-|beta-)?)(?P<mono>([dlDL]-)?(Glc|Gal|Man|Fuc|Xyl)([a-zA-Z5926,]+)?)(-)?(?P<link>(\([1-9?](/[1-9])*->([1-9?](/[1-9])*\))?)(-)?)(?P<bps>(\])*)"
    pattern = ""

    def patternCompilation(self):
        # Deprecated
        pa = r'(?P<anomer>(alpha-|beta-)?)'
        pm = r'(?P<mono>([dlDL]-)?(Glc|Gal|Man|Fuc|Xyl)([a-zA-Z5926,]+)?)(-)?'  # Why do I need 5926 here
        pl = r'(?P<link>(\([1-9?](/[1-9])*->([1-9?](/[1-9])*\))?)(-)?)'
        pc = r"(?P<bpe>(\[)*)" + pa + pm + pl + r"(?P<bps>(\])*)"
        compiledpattern = re.compile(pc)
        self.pattern = compiledpattern

    def regexSearch(self, seq):
        searchres = [m.groupdict() for m in self.pattern.finditer(seq)]
        searchres = list(reversed(searchres))

        res = []
        for s in searchres:
            res0 = copy.deepcopy(self.regexResTemplate)

            anomersl = s["anomer"]
            if anomersl == "alpha-":
                anomer = "Anomer.alpha"
            elif anomersl == "beta-":
                anomer = "Anomer.beta"
            elif anomersl == "":
                anomer = None
            else:
                # Unsupported anomer
                raise ValueError

            link = s["link"].lstrip("(").rstrip("-").rstrip(")")
            link = link.split("->")
            link = map(lambda x: None if x == "" else x, link)
            link = map(lambda x: None if x == "?" else x, link)

            res0["anomer"] = anomer
            res0["skel"] = s["mono"]
            res0["link"] = tuple(link)
            res0["substituent"] = None
            res0["sublink"] = None
            res0["branchpointstart"] = len(s["bps"])
            res0["branchpointend"] = len(s["bpe"])
        return searchres


class IUPACParserExtended1(IUPACParserAbstract):
    # Matthew's IUPAC, also the one in hayes column
    description = "UniCarb IUPAC string"
    example = "D-Neu5Ac(a2-3)[HSO3(-3)]D-Gal(b1-3)D-GalNAc(b1-4)[D-Neu5Ac(a2-8)D-Neu5Ac(a2-3)]D-Gal(b1-4)D-Glc(b1-1)Cer"
    alias = "IUPACExtended1:"
    precompiledpattern = r"(?P<bpe>(\[)*)(?P<sub>(\[[a-zA-Z0-9]+(\(-\d\))\]|[a-zA-Z0-9]+(\(-\d\))))?(?P<skel>([dlDL]-)?(Glc|Gal|Man|Fuc|Xyl|Neu|Ido)([a-zA-Z5926,]+)?)(?P<link>(\([ab]([1-9](/[1-9])*)-([1-9])?(/[1-9])*\)))?(?P<bps>(\])*)"

    def patternCompilation(self):
        # deprecated
        ps = r'(?P<sub>(\[[a-zA-Z0-9]+(\(-\d\))\]|[a-zA-Z0-9]+(\(-\d\))))?'
        pm = r'(?P<skel>([dlDL]-)?(Glc|Gal|Man|Fuc|Xyl|Neu|Ido|KDN)([a-zA-Z5926,]+)?)'
        pl = r'(?P<link>(\([ab][1-9]-([1-9])?\)))?'
        pattern = r"(?P<bpe>(\[)*)" + ps + pm + pl + r"(?P<bps>(\])*)"
        compiledpattern = re.compile(pattern)
        self.pattern = compiledpattern

    def regexSearch(self, seq):
        searchres = [m.groupdict() for m in self.pattern.finditer(seq)]
        searchres = list(reversed(searchres))
        if len(searchres) == 0:
            raise LinearCodeBadFormat(seq, 0)

        res = []
        for s in searchres:
            res0 = copy.deepcopy(self.regexResTemplate)

            if s["link"]:
                anomersl = s["link"][1]
                if anomersl == "a":
                    anomer = "Anomer.alpha"
                elif anomersl == "b":
                    anomer = "Anomer.beta"
                elif anomersl == "?":
                    anomer = None
                else:
                    raise IUPACLinearBadAnomer(seq, len(seq), anomersl)

                link = s["link"].lstrip("(").rstrip(")")
                link = link[1:].split("-")
                link = map(lambda x: None if x == "" else x, link)
                link = map(lambda x: None if x == "?" else x, link)
            else:
                anomer = None
                link = (None, None)

            skel = s["skel"]
            sub = s["sub"]
            # sublink = s["sublink"]
            if sub:
                sub = sub.lstrip("[").rstrip("]")
                skel = sub + s["skel"]
                # print sub

            res0["anomer"] = anomer
            res0["skel"] = skel
            res0["link"] = tuple(link)
            res0["substituent"] = sub
            # res0["sublink"] = sublink
            res0["branchpointstart"] = len(s["bps"])
            res0["branchpointend"] = len(s["bpe"])
            res.append(res0)
        return res


class IUPACParserCFG(IUPACParserAbstract):
    # The one in .cfg file
    description = ""
    example = "Fuca1-2Galb1-4(Fuca1-3)GlcNAcb1-3Galb1-4(Fuca1-3)GlcNAcb1-3Galb1-4(Fuca1-3)GlcNAcb-Sp0"
    alias = ""
    precompiledpattern = r"((?P<bpe>\()?)(?P<skel>((\([1-6SP]*\)|[1-6SP]*)*)?(Glc|Gal|Man|Fuc|Xyl|Neu|Ido|KDN|Rha|Mur)(([a-zA-Z5926,])*)?)(?P<link>\d?-?\d?)?(?P<bps>\))?"

    def regexSearch(self, seq):
        searchres = [m.groupdict() for m in self.pattern.finditer(seq)]
        searchres = list(reversed(searchres))
        if len(searchres) == 0:
            raise LinearCodeBadFormat(seq, 0)

        res = []
        for s in searchres:
            res0 = copy.deepcopy(self.regexResTemplate)

            skelori = self.skelreformat(s["skel"])
            skelwithoutanomer = s["skel"][:-1]
            anomer = s["skel"][-1]
            oriskelflag = False
            skelwithoutanomerflag = False
            anomerverification = anomer in "ab"
            try:
                self.mf.new(oriskelflag)
                oriskelflag = True
            except:
                pass

            try:
                self.mf.new(skelwithoutanomer)
                skelwithoutanomerflag = True
            except:
                pass

            if skelwithoutanomerflag and anomerverification:
                trueskel = skelwithoutanomer
                anomer = {"a": "Anomer.alpha", "b": "Anomer.beta"}[anomer]
            elif oriskelflag:
                trueskel = skelori
                anomer = None
            else:
                raise ValueError

            rawlink = s["link"]
            if rawlink:
                p1, p2 = rawlink.split("-")
                if not p1:
                    p1 = None
                if not p2:
                    p2 = None
                link = tuple([p1, p2])
            else:
                link = tuple([None, None])

            if not s["bps"]:
                s["bps"] = ""

            if not s["bpe"]:
                s["bpe"] = ""

            res0["anomer"] = anomer
            res0["skel"] = trueskel
            res0["link"] = tuple(link)
            res0["substituent"] = None
            res0["sublink"] = None
            res0["branchpointstart"] = len(s["bps"])
            res0["branchpointend"] = len(s["bpe"])
            res.append(res0)
        return res

    def skelreformat(self, skel):
        tri = r"(Glc|Gal|Man|Fuc|Xyl|Neu|Ido|KDN|Rha|Mur)"
        rest = re.compile(tri).split(skel)
        substr = rest.pop(0)
        rest = "".join(rest)

        reformattedsub = ""
        if substr:
            sub = re.compile(r"(\([1-6SP]+\)|[1-6SP]+)").findall(substr)
            sub = sorted(map(lambda x: x.lstrip("(").rstrip(")"), sub))

            if len(sub) == 1:
                reformattedsub = "(%s)" % sub[0]
            else:
                for i, s in enumerate(sub):
                    if i == 0:
                        reformattedsub += s
                    else:
                        reformattedsub += "(%s)" % s
        reformattedskel = reformattedsub + rest
        return reformattedskel


class IUPACParserGlycamExtended(IUPACParserAbstract):
    description = "The extended IUPAC format for Glycam.org"
    example = "DGalpa1-3DLyxpb1-5[DGalpNAcb1-3]DFrupb2-5DFrupb2-OME"
    alias = "Glycam:"
    precompiledpattern = r"(?P<bpe>(\[))?(?P<skel>((D|L)(Man|Gal|Glc|Ido|All|Alt|Gul|Tal|Xyl|Lyx|Rib|Ara|Fru|Psi|Sor|Tag|Fuc|Rha|Qui|KDN|KDO)(f|p)(NAc|A|5Ac|5Gc)?(a|b)))(?P<link>(\d-(\d)?))(?P<bps>(\]))?"
    
    # Not supported mono: Psi, Sor, Tag and Qui
    def regexSearch(self, seq):
        searchres = [m.groupdict() for m in self.pattern.finditer(seq)]
        searchres = list(reversed(searchres))
        if len(searchres) == 0:
            raise LinearCodeBadFormat(seq, 0)

        res = []
        for s in searchres:
            res0 = copy.deepcopy(self.regexResTemplate)

            skelwithoutanomer = s["skel"][:-1]
            anomer = s["skel"][-1]
            if anomer in "ab":
                anomer = {"a": "Anomer.alpha", "b": "Anomer.beta"}[anomer]
            else:
                raise ValueError

            rawlink = s["link"]
            if rawlink:
                p1, p2 = rawlink.split("-")
                if not p1:
                    p1 = None
                if not p2:
                    p2 = None
                link = tuple([p1, p2])
            else:
                link = tuple([None, None])

            if not s["bps"]:
                s["bps"] = ""

            if not s["bpe"]:
                s["bpe"] = ""

            res0["anomer"] = anomer
            res0["skel"] = skelwithoutanomer
            res0["link"] = tuple(link)
            res0["substituent"] = None
            res0["sublink"] = None
            res0["branchpointstart"] = len(s["bps"])
            res0["branchpointend"] = len(s["bpe"])
            res.append(res0)
        return res

class IUPACGlycamWriter:

    def __init__(self):
        self.sym = GlycamSym()

    def toString(self, g):
        nodes = list(g.all_nodes())
        nodesdict = {}
        relationship = {}
        linkinfo = {}
        for i, node in enumerate(nodes):
            nodesdict[i] = node

            if node.links():
                children = [x.child() for x in node.links()]
                children = map(lambda x: nodes.index(x), children)
            else:
                children = None
            relationship[i] = children

            if node.links():
                links = node.links()
                for l in links:
                    childindex = nodes.index(l.child())
                    pp = l._parent_pos
                    cp = l._child_pos
                    linkinfo[childindex] = tuple([pp, cp])

        # Debug purposes only
        #for k, v in enumerate(nodes):
        #    print k
        #    print v
        #    print

        monostrings = {}
        for i, node in nodesdict.items():
            if i == 0:
                monostrings[i] = self.mono2str(node, tuple([None, None]), root=True)
            else:
                monostrings[i] = self.mono2str(node, linkinfo[i], root=False)

        visited = ["[", "]"]
        tree = [0]
        flag = True
        while flag:
            flags = []
            for id in tree:
                if id in visited:
                    flags.append(True)
                    continue
                visited.append(id)
                children = relationship[id]
                if not children:
                    flags.append(True)
                    continue
                else:
                    formattedChildren = []
                    children = sorted(children, key=lambda x:self.ChildrenNumIncludingItself(x, nodesdict))
                    for i, c in enumerate(children):
                        if i == 0:
                            formattedChildren = [c] + formattedChildren
                        else:
                            formattedChildren = formattedChildren + ["[", c, "]"]
                    index = tree.index(id)
                    tree[index:index] = formattedChildren
                    flags.append(False)
                    break
            flag = False in flags
        newtree = []
        for monoid in tree:
            if monoid != "[" and monoid != "]":
                s = monostrings[monoid]
            else:
                s = monoid
            newtree.append(s)
        return "".join(newtree)

    def mono2str(self, m, l, root=False):
        sym = self.sym.toStr(m)

        config = "D"
        if m._config[0]:
            if m._config[0] == 2:
                config = "L"

        ringtype = "p"
        ringstart = m._ring_start
        ringend = m._ring_end
        if type(ringstart) == int or type(ringend) == int:
            if ringend - ringstart + 1 == 4:
                ringtype = "f"

        anomer = "?"
        if m._anomer:
            if m._anomer == 1:
                anomer = "a"
            elif m._anomer == 2:
                anomer = "b"

        skel = sym[0:3]
        modi = sym[3:]  # modifications, including mods and substituent

        if root:
            link = "1-"
        else:
            pp = l[0]
            cp = l[1]
            ppo, cpo = "", ""
            if pp == None:
                ppo = "?"
            else:
                pp = list(pp)
                pp = map(str, pp)
                ppo = "/".join(pp)

            if cp == None:
                cpo = "?"
            else:
                cp = list(cp)
                cp = map(str, cp)
                cpo = "/".join(cp)
            link = cpo + "-" + ppo
        return config + skel + ringtype + modi + anomer + link

    def ChildrenNumIncludingItself(self, id, d):
        gr = d[id]
        return -len(list(Glycan(gr).all_nodes()))

class IUPACGlycamFormat(GlycanFormatter):
    parser = IUPACParserGlycamExtended()
    writer = IUPACGlycamWriter()
    
    def toGlycan(self, seq):
        return self.parser.toGlycan(seq)
    
    def toStr(self, g):
        return self.writer.toString(g)

class WURCS20ParseError(GlycanParseError):
    pass


class WURCS20FormatError(WURCS20ParseError):
    def __init__(self,instr):
	self.message = "WURCS2.0 parser: Bad WURCS 2.0 format:\n      %s"%(instr,)

class UnsupportedMonoError(WURCS20ParseError):
    def __init__(self,monostr):
	self.message = "WURCS2.0 parser: Unsupported monosaccharide [%s]"%(monostr,)

class UnsupportedLinkError(WURCS20ParseError):
    def __init__(self,linkstr):
	self.message = "WURCS2.0 parser: Unsupported link %s"%(linkstr,)

class BadChildPositionLinkError(WURCS20ParseError):
    def __init__(self,linkstr):
	self.message = "WURCS2.0 parser: Bad child position in link %s"%(linkstr,)

class BadParentPositionLinkError(WURCS20ParseError):
    def __init__(self,linkstr):
	self.message = "WURCS2.0 parser: Bad parent position in link %s"%(linkstr,)

class MonoOrderLinkError(WURCS20ParseError):
    def __init__(self,linkstr):
	self.message = "WURCS2.0 parser: Unexpected monosaccharide order in link %s"%(linkstr,)

class CircularError(WURCS20ParseError):
  def __init__(self,instr):
   self.message = "WURCS2.0 parser: Circular structure not supported:\n     %s"%(instr,)

import WURCS20MonoFormatter

class WURCS20Format(GlycanFormatter):
    def __init__(self):
        self.mf = WURCS20MonoFormatter.WURCS20MonoFormat()
	self.wurcsre = re.compile(r'^WURCS=2\.0/(\d+,\d+,\d+)/((\[[^]]+\])+)/(\d+(-\d+)*)/(.*)$')
	self.simplelinkre = re.compile(r'^([a-zA-Z])([0-9?])-([a-zA-Z])([0-9?])$')
	self.multilinkre = re.compile(r'^([a-zA-Z])([0-9?])-(([a-zA-Z])([0-9?])(\|\4([0-9?]))*)$')
	self.ambiglinkre = re.compile(r'^([a-zA-Z])([0-9?])-(([a-zA-Z])([0-9?])(\|([a-zA-Z])([0-9?]))+)\}$')
        self.complinkre = re.compile(r'^([a-zA-Z][0-9?](\|[a-zA-Z][0-9?])*)\}-\{\1$')
        self.char2int = {}
        for i in range(26):
            self.char2int[chr(i+ord('a'))] = i+1
        for i in range(26):
            self.char2int[chr(i+ord('A'))] = i+26+1
    def toGlycan(self,s):
	m = self.wurcsre.search(s)
	if not m:
	    raise WURCS20FormatError(s)
	counts = map(int,m.group(1).split(','))
	distinctmono = {}; mono = {};
	for i,ms in enumerate(m.group(2)[1:-1].split('][')):
	    distinctmono[i+1] = ms
	for i,ms in enumerate(m.group(4).split('-')):
            mono[i+1] = self.mf.get(distinctmono[int(ms)])
            mono[i+1].set_id(i+1)

        undets = set()
        for li in map(str.strip,m.group(6).split('_')):

            if not li:
                continue
            
            mi = self.simplelinkre.search(li)
            if mi:
                
                ind1 = self.char2int[mi.group(1)]
                pos1 = (int(mi.group(2)) if mi.group(2) != "?" else None)
                ind2 = self.char2int[mi.group(3)]
                pos2 = (int(mi.group(4)) if mi.group(4) != "?" else None)

                if not (ind1 < ind2):
                    raise MonoOrderLinkError(li)

                if pos2 in (1,2,None):
                    parentmono = mono[ind1]
                    parentpos  = pos1
                    childmono  = mono[ind2]
                    childpos   = pos2
                elif pos1 in (1,) and pos2 not in (1,2,None):
                    # "Reversed" link?
                    parentmono = mono[ind2]
                    parentpos  = pos2
                    childmono  = mono[ind1]
                    childpos   = pos1

                # if not (childpos == None or childpos == childmono.ring_start()):
                #     raise BadChildPositionLinkError(li)

                # if not (parentpos == None or parentpos > parentmono.ring_start()):
                #     raise BadParentPositionLinkError(li)

                parentmono.add_child(childmono,
                                     child_pos=childpos,
                                     parent_pos=parentpos,
                                     parent_type=Linkage.oxygenPreserved,
                                     child_type=Linkage.oxygenLost)
                continue

            mi = self.multilinkre.search(li)
            if mi:

                ind1 = self.char2int[mi.group(1)]
                pos1 = (int(mi.group(2)) if mi.group(2) != "?" else None)
                ind2 = self.char2int[mi.group(4)]
                pos2 = sorted(set(map(lambda s: int(s[1]),mi.group(3).split('|'))))

                if not (ind2 < ind1):
                    raise MonoOrderLinkError(li)

                parentmono = mono[ind2]
                parentpos  = pos2
                childmono  = mono[ind1]
                childpos   = pos1

                # if not (childpos == None or childpos == childmono.ring_start()):
                #     raise BadChildPositionLinkError(li)

                # if not (parentpos[0] > parentmono.ring_start()):
                #     raise BadParentPositionLinkError(li)

                parentmono.add_child(childmono,
                                     child_pos=childpos,
                                     parent_pos=parentpos,
                                     parent_type=Linkage.oxygenPreserved,
                                     child_type=Linkage.oxygenLost)

                continue

            mi = self.ambiglinkre.search(li)
            if mi:

                ind1 = self.char2int[mi.group(1)]
                pos1 = (int(mi.group(2)) if mi.group(2) != "?" else None)

                indpos2 = defaultdict(set)
                for s in mi.group(3).split('|'):
                    ind2 = self.char2int[s[0]]
                    pos2 = (int(s[1]) if s[1] != '?' else None)
                    indpos2[ind2].add(pos2)
                for ind2 in indpos2:
                    indpos2[ind2] = sorted(indpos2[ind2])
                    if None in indpos2[ind2]:
                        indpos2[ind2] = None

		# print ind1,indpos2
                if not (max(indpos2) < ind1):
                    raise MonoOrderLinkError(li)

                if len(indpos2) < 2:
                    raise BadParentPositionLinkError(li)

                childmono  = mono[ind1]
                childpos   = pos1
                undets.add(childmono)

                # if not (childpos == None or childpos == 1 or childpos == childmono.ring_start()):
                #     raise BadChildPositionLinkError(li)

                for parentind,parentpos in indpos2.items():
                    parentmono = mono[parentind]
                    l = parentmono.add_child(childmono,
                                             child_pos=childpos,
                                             parent_pos=parentpos,
                                             parent_type=Linkage.oxygenPreserved,
                                             child_type=Linkage.oxygenLost)
                    l.set_undetermined(True)
                    l.set_instantiated(False)
		    l.child().set_connected(False)

                continue


            mi = self.complinkre.search(li)
            if mi:
                continue

            raise UnsupportedLinkError(li)

        unconnected = set()
        monocnt = 0
        for id,m in mono.items():
            monocnt += 1
            if len(m.parent_links()) > 0:
                continue
            m.set_connected(False)
            unconnected.add(m)
        # print len(unconnected)
        # for m in unconnected:
        #     print m
        if len(unconnected) == 0:
            raise CircularError(s)
        assert len(unconnected) in (1,monocnt), "# of unconnected nodes: %s not in {1,%d}"%(len(unconnected),monocnt)
        if len(unconnected) == 1:
            g = Glycan(unconnected.pop())
            g.set_undetermined(undets)
        else:
            assert len(undets) == 0
            g = Glycan()
            g.set_undetermined(unconnected)
	return g

if __name__ == '__main__':
    import sys, os.path
    clsinst = eval("%s()"%sys.argv[1])
    if len(sys.argv) > 2:
        files = sys.argv[2:]
    else:
        files = sys.stdin.read().splitlines()
    bad = 0
    for f in files:
        seq = open(f).read()
        try:
            g = clsinst.toGlycan(seq)
            print "+++", os.path.split(f)[1]
            # for t in g.undetermined_root_reprs():
            #     print t[1],str(t[0])
            print GlycoCTFormat().toStr(g)
        except GlycanParseError, e:
            print "!!!", os.path.split(f)[1], e
            bad += 1
    print "Failed: %d/%d"%(bad,len(files))
