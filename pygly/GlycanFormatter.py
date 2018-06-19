
from MonoFormatter import GlycoCTMonoFormat, LinCodeSym, LinCodeRank, IUPACSym
from Monosaccharide import Monosaccharide, Linkage, Anomer, Substituent, Mod
from Glycan import Glycan
from MonoFactory import MonoFactory
import re, sys
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

class GlycoCTUnexpectedLineError(GlycoCTParseError):
    def __init__(self,lineno,line):
	self.message = "GlycoCT parser, line %d: Unexpected line:\n ==> %s"%(lineno,line)
	self.line = line

class GlycoCTUndeterminedLinkageError(GlycoCTParseError):
    def __init__(self):
	self.message = "GlycoCT parser, undetermined linkage instantiation error"

class GlycoCTFormat(GlycanFormatter):
    def __init__(self):
        self.monofmt = GlycoCTMonoFormat()
    def toStr(self,g):
        g.set_ids()
        s = "RES\n"
        for m in g.subtree_nodes(g.root(),subst=True):
            s += self.monofmt.toStr(m)+"\n"
        s += "LIN\n"
        linkid = 1
        for l in sorted(g.subtree_links(g.root(),subst=True),key=lambda l: l.child().id()):
            l.set_id(linkid)
            linkid += 1
            s += self.monofmt.linkToStr(l)+"\n"
        undet = []
        for ur in g.undetermined_roots():
	    if not ur.connected():
	        undet.append(ur)
        if len(undet) == 0:
            return s
        s += "UND\n"
        for i,ur in enumerate(undet):
            s += "UND%s:100.0:100.0\n"%(i+1,)
            parentids = set()
            subtreelinkage = set()
            for pl in ur.parent_links():
                parentids.add(pl.parent().id())
                subtreelinkage.add(self.monofmt.linkToStr(pl,noids=True))
            s += "ParentIDs:%s\n"%("|".join(map(str,sorted(parentids))),)
            assert len(subtreelinkage) == 1
            s += "SubtreeLinkageID1:%s\n"%(subtreelinkage.pop(),)
            s += "RES\n"
            for m in g.subtree_nodes(ur,subst=True):
                s += self.monofmt.toStr(m)+"\n"
            first = True
            for l in sorted(g.subtree_links(ur,subst=True),key=lambda l: l.child().id()):
                if first:
                    s += "LIN\n"
                    first = False
                l.set_id(linkid)
                linkid += 1
                s += self.monofmt.linkToStr(l)+"\n"
        return s

    def toGlycan(self,s):
        res = {}
	und = defaultdict(dict)
	undind = None
	state = None
	undet = False
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
			undet = True
			for l in links:
			    l.set_undetermined(True)
			    l.set_instantiated(False)
			    l.child().set_connected(False)
		except RuntimeError, e:
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
	    undet = True
	    linkagestr = d['stlink'][1]
	    rootid = d['root']
	    for pid in d['parentids']:
		links = self.monofmt.linkFromStr("0:%s%s%s%s"%(pid,linkagestr[:-1],rootid,linkagestr[-1]),res)
		for l in links:
		    l.set_undetermined(True)
		    l.set_instantiated(False)
		    l.child().set_connected(False)
	g = Glycan(res[1])
	g.set_undetermined(undets)
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

class WURCS20ParseError(GlycanParseError):
    pass


class WURCS20FormatError(WURCS20ParseError):
    def __init__(self,instr):
	self.message = "WURCS2.0 parser: Bad WURCS 2.0 format %s"%(instr,)

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

class WURCS20Format(GlycanFormatter):
    def __init__(self):
        self.mf = MonoFactory()
	self.wurcsre = re.compile(r'^WURCS=2\.0/(\d+,\d+,\d+)/((\[[^]]+\])+)/(\d+(-\d+)*)/(.*)$')
	self.simplelinkre = re.compile(r'^([a-zA-Z])([0-9?])-([a-zA-Z])([0-9?])$')
	self.multilinkre = re.compile(r'^([a-zA-Z])([0-9?])-(([a-zA-Z])([0-9?])(\|\4([0-9?]))*)$')
	self.ambiglinkre = re.compile(r'^([a-zA-Z])([0-9?])-(([a-zA-Z])([0-9?])(\|([a-zA-Z])([0-9?]))+)\}$')
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
            try:
                mono[i+1] = self.mf.new('WURCS20:['+distinctmono[int(ms)]+']')
            except KeyError:
                raise UnsupportedMonoError(distinctmono[int(ms)])
	root = mono[1]
        undet = False
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

                parentmono = mono[ind1]
                parentpos  = pos1
                childmono  = mono[ind2]
                childpos   = pos2

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

                if len(parentpos) == 2:
                    parentmono.add_child(childmono,
                                         child_pos=childpos,
                                         parent_pos=parentpos[0],
                                         parent_pos2=parentpos[1],
                                         parent_type=Linkage.oxygenPreserved,
                                         child_type=Linkage.oxygenLost)
                else:
                    BadParentPositionLinkError(li)

                continue

            mi = self.ambiglinkre.search(li)
            if mi:

                undet = True

                ind1 = self.char2int[mi.group(1)]
                pos1 = (int(mi.group(2)) if mi.group(2) != "?" else None)

                indpos2 = defaultdict(set)
                for s in mi.group(3).split('|'):
                    ind2 = self.char2int[s[0]]
                    pos2 = (int(s[1]) if s[1] != '?' else None)
                    indpos2[ind2].add(pos2)
                for ind2 in indpos2:
                    indpos2[ind2] = sorted(indpos2[ind2])

		# print ind1,indpos2
                if not (max(indpos2) < ind1):
                    raise MonoOrderLinkError(li)

                if len(indpos2) < 2:
                    raise BadParentPositionLinkError(li)

                for ind2 in indpos2:
                    if len(indpos2[ind2]) > 2:
                        raise BadParentPositionLinkError(li)
                    if None in indpos2[ind2] and len(indpos2[ind2]) != 1:
                        raise BadParentPositionLinkError(li)

                childmono  = mono[ind1]
                childpos   = pos1
                undets.add(childmono)

                # if not (childpos == None or childpos == 1 or childpos == childmono.ring_start()):
                #     raise BadChildPositionLinkError(li)

                for parentind,parentpos in indpos2.items():
                    parentmono = mono[parentind]
                    parentpos2 = None
                    if len(parentpos) == 2:
                        parentpos2=parentpos[1]
                    parentpos=parentpos[0]
                    l = parentmono.add_child(childmono,
                                             child_pos=childpos,
                                             parent_pos=parentpos,
                                             parent_pos2=parentpos2,
                                             parent_type=Linkage.oxygenPreserved,
                                             child_type=Linkage.oxygenLost)
                    l.set_undetermined(True)
                    l.set_instantiated(False)
		    l.child().set_connected(False)
                continue

##             mi = self.ambiglinkre.search(li)
##             if mi:
##                 undet = True
##                 tomono = mono[self.char2int[mi.group(1)]]
##                 if mi.group(2) == "?":
##                     inpos = None
##                 else:
##                     inpos = int(mi.group(2))
##                 for plink in mi.group(3).split('|'):

            raise UnsupportedLinkError(li)

        g = Glycan(mono[1])
        g.set_undetermined(undets)
        g.set_ids()
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
        except GlycanParseError, e:
            print "!!!", os.path.split(f)[1], e
            bad += 1
    print "Failed: %d/%d"%(bad,len(files))
