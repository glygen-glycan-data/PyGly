
from __future__ import print_function

from MonoFormatter import GlycoCTMonoFormat, LinCodeSym, LinCodeRank, IUPACSym, GlycamSym
from Monosaccharide import Monosaccharide, Linkage, Anomer, Substituent, Mod
from Glycan import Glycan
from MonoFactory import MonoFactory

import re, sys, traceback
import copy
import string
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

class GlycoCTREPBeforeRESLINError(GlycoCTParseError):
    def __init__(self,lineno):
        self.message = "GlycoCT parser, line %d: REP section before RES or LIN section"%(lineno,)

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

class GlycoCTRepeatTimesError(GlycoCTParseError):
    def __init__(self):
        self.message = "GlycoCT parser, repeat times range max is smaller than min"

class GlycoCTRepeatInRepeatError(GlycoCTParseError):
    def __init__(self):
        self.message = "GlycoCT parser, found repeat in repeat"

class GlycoCTUnconnectedCountError(GlycoCTParseError):
    def __init__(self):
        self.message = "GlycoCT parser, unconnected node count error"

class GlycoCTFormat(GlycanFormatter):
    def __init__(self):
        self.monofmt = GlycoCTMonoFormat()
    
    def monoteardown(self, m):
        m = copy.deepcopy(m)
        m.set_id(None)
        mstr = self.monofmt.toStr(m)
        mstr = mstr[2:]

        subs = []
        for sl in list(m.substituent_links()):

            if sl.parent_pos():
                parent_pos = sorted(list(sl.parent_pos()))
                pp = "/".join(map(str, parent_pos))
            else:
                pp = "?"

            if sl.child_pos():
                child_pos = sorted(list(sl.child_pos()))
                cp = "/".join(map(str, child_pos))
            else:
                cp = "?"

            if sl.parent_type():
                pt = {1: "o",
                      2: "d",
                      3: "x",
                      4: "n"}.get(list(sl.parent_type())[0], "?")
            else:
                pt = "?"

            sl.child().set_id(None)
            subname = self.monofmt.toStr(sl.child())[2:]

            subs.append({"pp": pp,
                         "cp": cp,
                         "pt": pt,
                         "subname": subname, })

        return (mstr, subs)

    def mtoStr(self, m):
        mstr, subs = self.monoteardown(m)

        substrs = []
        for s in subs:
            pp = s["pp"]
            cp = s["cp"]
            pt = s["pt"]
            subname = s["subname"]

            substr = "(%s%s:%s)%s" % (pp, pt, cp, subname)
            substrs.append(substr)

        if substrs:
            substrs = sorted(substrs)
            allsub = "|".join(substrs)
            mstr += "||%s" % allsub

        return mstr

    def mtodict(self, m):
        mstr, subs = self.monoteardown(m)

        subdict = []
        for s in subs:
            pp = s["pp"]
            cp = s["cp"]
            pt = s["pt"]
            subname = s["subname"]

            d = {"substMsPos": pp,
                 "substName": subname,
                 "substMsLinktype": pt}
            substr = "(%s%s:%s)%s" % (pp, pt, cp, subname)
            subdict.append(d)

        return (mstr, subdict)
    
    def toStr(self, g):

        g.unset_ids()


        all_nodes = list(g.all_nodes(subst=True))

        # Modify the entry and exit link before output the string
        entry_link_child_type_backup = "Default value"
        exit_link_parent_type_backup = "Default value"
        for node in all_nodes:

            for l in node.links():

                # Entry link
                if l.child().is_repeat_start():
                    entry_link_child_type_backup = l._child_type
                    l.set_child_type(Linkage.nitrogenAdded)

                # Exit link
                if l.repeat_unit_out_link():
                    exit_link_parent_type_backup = l._parent_type
                    l.set_parent_type(Linkage.nitrogenAdded)


        repeat_info = {}
        glycoct_node_id_to_rep_id = {}
        rep_id_to_glycoct_node_id = {}



        # Temporary stuff
        repeat_starts = []
        repeat_starts_d = {}
        nodes_in_repeat = []
        for node in all_nodes:
            if node.is_repeat_start():
                repeat_starts.append(node)


        for i, rs in enumerate(repeat_starts):
            rep_index = i+1

            repeat_starts_d[rs] = rep_index
            for pl in rs.parent_links(include_repeat=True):
                if pl.repeat_bridge_link():
                    # repeat_start, repeat_end, nodeset, linkset
                    repeat_info[rep_index] = [None, None, [], []]

                    re = pl.parent()
                    repeat_info[rep_index][0] = rs
                    repeat_info[rep_index][1] = re


                    todo = [rs]

                    while len(todo) > 0:
                        n = todo.pop()

                        if n not in repeat_info[rep_index][2]:
                            repeat_info[rep_index][2].append(n)


                        substs = n.substituents()
                        todo += substs

                        if n == re:
                            inside_repeat_unit_child = []
                            for l in re.links():
                                if l.basic_link():
                                    inside_repeat_unit_child.append(l.child())
                            todo += inside_repeat_unit_child
                        else:
                            todo += n.children()

                        #repeat_info[rep_index][2]

                    nodes_in_repeat += repeat_info[rep_index][2]
                    # print(len(repeat_info[rep_index][2]))


            assert None not in repeat_info[rep_index]



        i = 1
        for n in all_nodes:

            if n in nodes_in_repeat:

                if n in repeat_starts_d:
                    glycoct_node_id_to_rep_id[i] = repeat_starts_d[n]
                    rep_id_to_glycoct_node_id[repeat_starts_d[n]] = i
                    i += 1
                continue

            n.set_id(i)
            i += 1

        for n in all_nodes:

            if n.id() == None:
                n.set_id(i)
                i += 1


        roots = []
        if g.root():
            roots.append(g.root())
        for ur in g.undetermined_roots():
            if not ur.connected():
                roots.append(ur)



        s = "RES\n"


        for r in roots:
            if len(r.parent_links()) == 0 or r == g.root():
                for m in g.subtree_nodes(r, subst=True):

                    if m in repeat_starts_d:
                        rep_id = repeat_starts_d[m]
                        glycoct_node_id = rep_id_to_glycoct_node_id[rep_id]
                        s += "%sr:r%s\n" % (glycoct_node_id, rep_id)

                    if m in nodes_in_repeat:
                        continue

                    s += self.monofmt.toStr(m).strip('!')+"\n"


        first = True
        linkid = 1
        for r in roots:
            if len(r.parent_links()) == 0 or r == g.root():

                sortkey = {}
                for rep_id, rinfo in repeat_info.items():
                    glycoct_node_id = rep_id_to_glycoct_node_id[rep_id]
                    rs = rinfo[0]
                    sortkey[rs.id()] = glycoct_node_id

                for l in sorted(g.subtree_links(r,subst=True),key=lambda l: l.child().id() if l.child().id() not in sortkey else sortkey[l.child().id()] ):

                    if l.parent() in nodes_in_repeat or l.child() in nodes_in_repeat:

                        lstr = None
                        if l.parent() not in nodes_in_repeat and l.child() in nodes_in_repeat:
                            # Entry to repeat unit
                            l.set_id(linkid)
                            linkid += 1

                            rep_id_found = None
                            for rep_id, rinfo in repeat_info.items():
                                if l.child() == rinfo[0]:
                                    rep_id_found = rep_id
                            assert rep_id_found != None
                            replace_child_id = rep_id_to_glycoct_node_id[rep_id_found]

                            lstr = self.monofmt.linkToStr(l, replace_child_id=replace_child_id)


                        if l.parent() in nodes_in_repeat and l.child() not in nodes_in_repeat:
                            # Exit repeat unit
                            l.set_id(linkid)
                            linkid += 1

                            rep_id_found = None
                            for rep_id, rinfo in repeat_info.items():
                                if l.parent() == rinfo[1]:
                                    rep_id_found = rep_id
                            assert rep_id_found != None
                            replace_parent_id = rep_id_to_glycoct_node_id[rep_id_found]

                            lstr = self.monofmt.linkToStr(l, replace_parent_id=replace_parent_id)


                        if lstr != None:
                            if first:
                                s += "LIN\n"
                                first = False
                            s += lstr + "\n"
                        else:
                            if l.repeat_bridge_link():
                                continue

                            found = False
                            for rep_id, rinfo in repeat_info.items():

                                if l.child() in rinfo[2] and l.parent() in rinfo[2]:
                                    rinfo[3].append(l)
                                    found = True
                            assert found

                        continue

                    if first:
                        s += "LIN\n"
                        first = False
                    l.set_id(linkid)
                    linkid += 1
                    s += self.monofmt.linkToStr(l)+"\n"


        if len(repeat_info) > 0:
            s += "REP\n"

            for rep_id in sorted(repeat_info.keys()):
                rinfo = repeat_info[rep_id]
                rs = rinfo[0]

                repeat_bridge_link = filter(lambda pl: pl.repeat_bridge_link(), rs.parent_links(include_repeat=True))[0]
                rmin = repeat_bridge_link.repeat_times_min()
                if rmin == None:
                    rmin = -1
                rmax = repeat_bridge_link.repeat_times_max()
                if rmax == None:
                    rmax = -1

                s += "REP%s:%s=%s-%s\n" % (rep_id, self.monofmt.linkToStr(repeat_bridge_link), rmin, rmax)

                s += "RES\n"
                for n in sorted(rinfo[2], key=lambda n: n.id()):
                    s += self.monofmt.toStr(n).strip('!') + "\n"

                if len(rinfo[3]) > 0:
                    s += "LIN\n"

                    for l in rinfo[3]:
                        l.set_id(linkid)
                        linkid += 1
                        s += self.monofmt.linkToStr(l) + "\n"


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


        for node in all_nodes:

            for l in node.links():

                # Entry link
                if l.child().is_repeat_start():
                    assert entry_link_child_type_backup != "Default value"
                    l.set_child_type(entry_link_child_type_backup)

                # Exit link
                if l.repeat_unit_out_link():
                    assert exit_link_parent_type_backup != "Default value"
                    l.set_parent_type(exit_link_parent_type_backup)

        return s

    def toGlycan(self,s):
        res = {}

        rep = {}

        und = defaultdict(dict)
        undind = None
        undets = set()

        state = ""
        seen = set()

        for lineno,l in enumerate(s.splitlines()):
            l = l.strip()
            if not l:
                continue

            if l == "RES":
                if l in seen and ("UND" not in state and "REP" not in state):
                    raise GlycoCTRepeatedSectionError(lineno=lineno+1,section=l)
                if state == "UND":
                    state = "UNDRES"
                elif "REP" in state:
                    state = "REPRES"
                else:
                    state = "RES"
                seen.add(state)
                continue

            if l == "LIN":
                if l in seen and ("UND" not in state and "REP" not in state):
                    raise GlycoCTRepeatedSectionError(lineno=lineno+1,section=l)
                if "RES" not in seen:
                    raise GlycoCTLINBeforeRESError()
                if "UND" in state:
                    state = "UNDLIN"
                elif "REP" in state:
                    state = "REPLIN"
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

            if l == "REP":
                if l in seen:
                    raise GlycoCTRepeatedSectionError(lineno=lineno+1,section=l)
                if "RES" not in seen:
                    raise GlycoCTREPBeforeRESLINError()
                state = "REP"
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

            if state in ("RES", "UNDRES", "REPRES"):

                x = re.findall(r'^(\d+)r:r(\d+)$', l)
                if len(x) == 1:
                    glycoctid, glycoctrepeatid = map(int, x[0])
                    rep[glycoctid] = {"id": glycoctrepeatid}
                    continue

                try:
                    m = self.monofmt.fromStr(l)
                except (TypeError,LookupError,ValueError,RuntimeError) as e:
                    raise GycoCTBadRESLineError(message=e.args[0],lineno=lineno+1,line=l)
                res[m.id()] = m
                if state == "UNDRES" and undind != None and 'root' not in und[undind]:
                    und[undind]['root'] = m.id()
                    undets.add(m)
                continue

            if state in ("LIN", "REP", "REPLIN", "UNDLIN"):
                id, parentid, parentpos, parenttype, childid, childpos, childtype, r1, r2, r3 = self.monofmt.linkParseBase(l)

                if parentid in rep:
                    rep[parentid]["link_out"] = l
                    continue

                if childid in rep:
                    rep[childid]["link_in"] = l
                    continue

                if r1:
                    # Repeat bridge
                    tmp = None
                    for glycoctid in rep.keys():
                        glycoctrepeatid = rep[glycoctid]["id"]
                        if glycoctrepeatid == id:
                            tmp = glycoctid
                    assert tmp != None

                    rep[tmp]["start"] = childid
                    rep[tmp]["end"] = parentid
                    rep[tmp]["repeat_bridge"] = l

                    rmin = 0
                    if r2 != None:
                        rmin = r2

                    rmax = 2e100
                    if r3 != None:
                        rmax = r3

                    if rmin > rmax:
                        raise GlycoCTRepeatTimesError

                    rep[tmp]["repeat_times"] = (r2, r3)
                    continue

                try:
                    links = self.monofmt.linkFromStr(l, res)
                    if len(links) > 1:
                        for l in links:
                            l.set_undetermined(True)
                            l.set_instantiated(False)
                            l.child().set_connected(False)
                except (RuntimeError,AttributeError) as e:
                    raise GlycoCTBadLINLineError(message=e.args[0],lineno=lineno+1,line=l)      
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
                except (RuntimeError,AttributeError) as e:
                    raise GlycoCTParentLinkError(message=e.args[0],lineno=lineno+1,line=l)


        unconnected = set()

        for glycoctid in rep:

            assert isinstance(res[rep[glycoctid]["start"]] ,Monosaccharide)
            assert isinstance(res[rep[glycoctid]["end"]], Monosaccharide)

            entry_link_child_type = None
            exit_link_parent_type = None
            assert "repeat_bridge" in rep[glycoctid]
            if "repeat_bridge" in rep[glycoctid]:
                id, parentid, parentpos, parenttype, childid, childpos, childtype, r1, r2, r3 = self.monofmt.linkParseBase(rep[glycoctid]["repeat_bridge"])
                parentid = rep[glycoctid]["end"]
                childid = rep[glycoctid]["start"]
                parent = res[parentid]
                child = res[childid]

                entry_link_child_type = childtype
                exit_link_parent_type = parenttype

                l = self.monofmt.linkFromPara(parent, child, parenttype, parentpos, childtype, childpos, repeat_bridge=True)

                l.set_repeat_times_min(rep[glycoctid]["repeat_times"][0])
                l.set_repeat_times_max(rep[glycoctid]["repeat_times"][1])

            if "link_in" in rep[glycoctid]:
                id, parentid, parentpos, parenttype, childid, childpos, childtype, r1, r2, r3 = self.monofmt.linkParseBase(rep[glycoctid]["link_in"])
                childid = rep[glycoctid]["start"]
                parent = res[parentid]
                child = res[childid]
                self.monofmt.linkFromPara(parent, child, parenttype, parentpos, entry_link_child_type, childpos)
            else:
                r = res[rep[glycoctid]["start"]]
                unconnected.add(r)

            if "link_out" in rep[glycoctid]:
                id, parentid, parentpos, parenttype, childid, childpos, childtype, r1, r2, r3 = self.monofmt.linkParseBase(rep[glycoctid]["link_out"])
                parentid = rep[glycoctid]["end"]
                parent = res[parentid]
                child = res[childid]
                self.monofmt.linkFromPara(parent, child, exit_link_parent_type, parentpos, childtype, childpos, repeat_exit=True)


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

        if g.is_repeat_in_repeat():
            raise GlycoCTRepeatInRepeatError()

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


class IUPACParseError(GlycanParseError):
    pass

class IUPACBranchingError(IUPACParseError):
    def __init__(self):
        self.message = "Input IUPAC string contains multiple adjacent branch point, which is not allowed"

class IUPACUnsupportedAnomer(IUPACParseError):
    def __init__(self, sym, a):
        self.message = "Unsupported anomer: %s in link/mono: %s" % (a, sym)

class IUPACUnsupportedSym(IUPACParseError):
    def __init__(self, sym):
        self.message = "Unsupported sym: %s" % sym

class IUPACSkippedMonosaccharide(IUPACParseError):
    def __init__(self):
        self.message = "Parsing error, some monosaccharide from IUPAC string are not parsed in"


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
                        raise IUPACBranchingError
                    branchPoint.append(parentMono)
                parentMono = mono

                # Well... It can make sense if there are multiple ending bracket. Just pop n element from branchPoint
                if remove_stack != 0:
                    if remove_stack > 1:
                        raise IUPACBranchingError
                    parentMono = branchPoint.pop(-1)
        glycanObj = Glycan(root)
        if len(list(glycanObj.all_nodes())) != self.monoNumCheck(seq):
            raise IUPACSkippedMonosaccharide
        return glycanObj

    def monoNumCheck(self, seq):
        possibleTriletter = "Man|Gal|Glc|Ido|All|Alt|Gul|Tal|Xyl|Lyx|Rib|Ara|Fru|Psi|Sor|Tag|Fuc|Rha|Qui|KDN|KDO|Neu"
        possibleTriletter = possibleTriletter.lower()
        ptl = possibleTriletter.split("|")
        s = seq.lower()
        c = 0
        for tri in ptl:
            c += s.count(tri)
        return c

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
                raise IUPACUnsupportedAnomer(s["link"], anomersl)

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
                raise IUPACUnsupportedAnomer(s["mono"], anomersl)

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
    precompiledpattern = r"(?P<bpe>(\[)*)(?P<sub>(\[[a-zA-Z0-9]+(\(-\d\))\]|[a-zA-Z0-9]+(\(-\d\))))?(?P<skel>([dlDL]-)?(Glc|Gal|Man|Fuc|Xyl|Neu|Ido|KDN)([a-zA-Z5926,]+)?)(?P<link>(\([ab]([1-9](/[1-9])*)-([1-9])?(/[1-9])*\)))?(?P<bps>(\])*)"

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
                    raise IUPACUnsupportedAnomer(s["link"], anomersl)

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
                raise IUPACUnsupportedSym(skelori)

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
    precompiledpattern = r"(?P<bpe>(\[))?(?P<skel>((D|L)(Man|Gal|Glc|Ido|All|Alt|Gul|Tal|Xyl|Lyx|Rib|Ara|Fru|Psi|Sor|Tag|Fuc|Rha|Qui|KDN|KDO)(f|p)(NAc|A|5Ac|5Gc)?(\[(((\d(S|P|A|Me))(,)?)+)\])?(a|b)))(?P<link>(\d-(\d)?))(?P<bps>(\]))?"
    
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
                raise IUPACUnsupportedAnomer(s["skel"], ["skel"][-1])

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

    def has1v1children(self, m):

        for cl in m.links():
            cp = cl._child_pos
            pp = cl._parent_pos
            if cp == {1} and pp == {1}:
                return True
        return False

    def mono2str(self, m, l, root=False):
        m, subsstr = self.subs2str(m)
        sym = self.sym.toStr(m)

        config = "?"
        if m._config:
            if m._config[0] == 1:
                config = "D"
            if m._config[0] == 2:
                config = "L"

        ringtype = "?"
        ringstart = m._ring_start
        ringend = m._ring_end
        if type(ringstart) == int and type(ringend) == int:
            if ringend - ringstart + 1 == 4:
                ringtype = "f"
            elif ringend - ringstart + 1 == 5:
                ringtype = "p"

        anomer = "?"
        if m._anomer:
            if m._anomer == 1:
                anomer = "a"
            elif m._anomer == 2:
                anomer = "b"

        skel = sym[0:3]
        modi = sym[3:]  # modifications, including mods and substituent

        if root:
            link = "1-OH"

            # Glycam has more strict layout for 1v1 case
            assert not self.has1v1children(m)
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
        return config + skel + ringtype + modi + subsstr + anomer + link

    def subs2str(self, m):
        m = copy.deepcopy(m)
        subs = []
        sub2str = {
            Substituent.sulfate: "S",
            Substituent.phosphate: "P",
            Substituent.acetyl: "A",
            Substituent.methyl: "Me"
        }

        for sub in m.substituent_links():
            if sub._parent_pos == None:
                pp = "?"
            else:
                pp = "|".join(map(str, list(sub._parent_pos)))

            if sub.child()._sub != Substituent.nAcetyl:
                try:
                    subname = sub2str[sub.child()._sub]
                except KeyError:
                    raise KeyError
                subs.append(pp + subname)

                m._substituent_links.remove(sub)

        subs = sorted(subs)
        if subs:
            subsstr = "[" + "/".join(subs) + "]"
        else:
            subsstr = ""
        return m, subsstr

    def ChildrenNumIncludingItself(self, id, d):
        gr = d[id]
        num = len(list(Glycan(gr).all_nodes()))
        if num == 1 and gr._stem == tuple([14]) and gr._mods == [((6,), 1)]:
            # this method is used to sort the children branches
            # If we have multiple branch under a branch point
            # And one branch is a single Fructose
            # Then it makes more sense to put this single fructose in a "[]"
            return 0
        return -num

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

class UndeterminedLinkCountError(WURCS20ParseError):
    def __init__(self):
        self.message = "WURCS2.0 parser: undetermined link count"

class ZeroPlusLinkCountError(WURCS20ParseError):
    def __init__(self):
        self.message = "WURCS2.0 parser: zero plus link count"

class BadChildPositionLinkError(WURCS20ParseError):
    def __init__(self,linkstr):
        self.message = "WURCS2.0 parser: Bad child position in link %s"%(linkstr,)

class BadParentPositionLinkError(WURCS20ParseError):
    def __init__(self,linkstr):
        self.message = "WURCS2.0 parser: Bad parent position in link %s"%(linkstr,)

class BadRepeatTimesError(WURCS20ParseError):
    def __init__(self,linkstr):
        self.message = "WURCS2.0 parser: Bad repeat times range in link %s (max is smaller than min)" % (linkstr,)

class MonoOrderLinkError(WURCS20ParseError):
    def __init__(self,linkstr):
        self.message = "WURCS2.0 parser: Unexpected monosaccharide order in link %s"%(linkstr,)

class LinkCountError(WURCS20ParseError):
  def __init__(self,instr):
   self.message = "WURCS2.0 parser: Bad link count:\n     %s"%(instr,)

class CircularError(WURCS20ParseError):
  def __init__(self,instr):
   self.message = "WURCS2.0 parser: Circular structure not supported:\n     %s"%(instr,)

class UnexpectedConnectivityError(WURCS20ParseError):
  def __init__(self,instr):
   self.message = "WURCS2.0 parser: Unexpected or strange connectivity:\n     %s"%(instr,)

class UnexpectedFloatingSubstError(WURCS20ParseError):
  def __init__(self,instr):
   self.message = "WURCS2.0 parser: Unexpected floating substituent:\n     %s"%(instr,)

class UnexpectedRepeatInRepeatError(WURCS20ParseError):
  def __init__(self):
   self.message = "WURCS2.0 parser: Unexpected found (repeat in repeat)"

import WURCS20MonoFormatter

class WURCS20Format(GlycanFormatter):
    def __init__(self):
        self.mf = WURCS20MonoFormatter.WURCS20MonoFormat()
        self.wurcsre = re.compile(r'^WURCS=2\.0/(\d+,\d+,\d+\+?)/((\[[^]]+\])+)/(\d+(-\d+)*)/(.*)$')
        self.simplelinkre = re.compile(r'^([a-zA-Z]{1,2})([0-9?])-([a-zA-Z]{1,2})([0-9?])$')
        self.multilinkre = re.compile(r'^([a-zA-Z]{1,2})([0-9?])-(([a-zA-Z]{1,2})([0-9?])(\|\4([0-9?]))*)$')
        self.ambiglinkre = re.compile(r'^([a-zA-Z]{1,2})([0-9?])-(([a-zA-Z]{1,2})([0-9?])(\|([a-zA-Z]{1,2})([0-9?]))+)\}$')
        self.complinkre = re.compile(r'^([a-zA-Z]{1,2}[0-9?](\|[a-zA-Z]{1,2}[0-9?])*)\}-\{\1$')
        self.compsubstlinkre = re.compile(r'^([a-zA-Z]{1,2}[0-9?](\|[a-zA-Z]{1,2}[0-9?])*)\}\*(.*)$')
        self.substbridgelinkre = re.compile(r"^([a-zA-Z]{1,2})([0-9?])-([a-zA-Z]{1,2})([0-9?])\*(.*?)$")
        self.substbridgecompre = re.compile(r"^([a-zA-Z]{1,2}[?](\|[a-zA-Z]{1,2}[?])*)\}-\{([a-zA-Z]{1,2}[?](\|[a-zA-Z]{1,2}[?])*)\*(.*?)$")
        self.repeatbridgelinkre = re.compile(r'^([a-zA-Z]{1,2})([0-9?])-([a-zA-Z]{1,2})([0-9?])(\*.*)?~(n|\d+:(\d+|n))$')

        self.char2int = {}
        self.int2char = {}
        for i in range(26):
            self.char2int[chr(i+ord('a'))] = i+1
            self.int2char[i+1] = chr(i+ord('a'))
        for i in range(26):
            self.char2int[chr(i+ord('A'))] = i+26+1
            self.int2char[i + 1 + 26] = chr(i+ord('A'))

        allLetters = string.ascii_lowercase + string.ascii_uppercase
        for l1 in allLetters:
            for l2 in allLetters:
                doubledigit = l1+l2
                intfordoubledigit = self.char2int[l2] + self.char2int[l1]*52
                self.int2char[intfordoubledigit] = doubledigit
                self.char2int[doubledigit] = intfordoubledigit

    def toGlycan(self,s):
        m = self.wurcsre.search(s)
        if not m:
            raise WURCS20FormatError(s)
        if m.group(1).endswith('+'):
            if m.group(1).endswith(",0+"):
                raise ZeroPlusLinkCountError()
            else:
                raise UndeterminedLinkCountError()

        counts = map(int,m.group(1).split(','))
        distinctmono = {}
        mono = {}
        rep = []

        for i,ms in enumerate(m.group(2)[1:-1].split('][')):
            distinctmono[i+1] = ms
        for i,ms in enumerate(m.group(4).split('-')):
            mono[i+1] = self.mf.get(distinctmono[int(ms)])
            mono[i+1].set_id(i+1)
            mono[i+1].set_external_descriptor_id(i+1)

        undets = set()
        floating_substs = []

        # Note, repeat linkage is always at the end of link_list
        # Need some info from repeat bridge first
        wurcs_linkage_list = [ s.strip() for s in m.group(6).split('_') ]
        wurcs_repeat_linkage_list = []
        for l in reversed(wurcs_linkage_list):
            if "~" in l:
                wurcs_linkage_list.remove(l)
                wurcs_repeat_linkage_list.append(l)


        for li in wurcs_repeat_linkage_list + wurcs_linkage_list:

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
                elif pos1 in (1,None) and pos2 not in (1,2,None):
                    # "Reversed" link?
                    parentmono = mono[ind2]
                    parentpos  = pos2
                    childmono  = mono[ind1]
                    childpos   = pos1
                else:
                    raise UnsupportedLinkError(li)

                # if not (childpos == None or childpos == childmono.ring_start()):
                #     raise BadChildPositionLinkError(li)

                # if not (parentpos == None or parentpos > parentmono.ring_start()):
                #     raise BadParentPositionLinkError(li)

                self.add_child(parentmono, childmono, parentpos, childpos, rep)
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

                self.add_child(parentmono, childmono, parentpos, childpos, rep)
                continue

            mi = self.substbridgecompre.search(li)
            if mi:
                # composition like substituent in link-situation

                subst = self.mf.getsubst(mi.group(5))
                subst.set_connected(False)
                floating_substs.append(subst)
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

            mi = self.compsubstlinkre.search(li)
            if mi:
                subst = self.mf.getsubst(mi.group(3))
                subst.set_connected(False)
                floating_substs.append(subst)
                continue


            mi = self.repeatbridgelinkre.search(li)
            if mi:
                start_ind = self.char2int[mi.group(1)]
                start_pos = (int(mi.group(2)) if mi.group(2) != "?" else None)
                end_ind = self.char2int[mi.group(3)]
                end_pos = (int(mi.group(4)) if mi.group(4) != "?" else None)

                wurcssubststr = mi.group(5)

                repeattimes_raw = mi.group(6)
                repeattimes = (None, None)
                if repeattimes_raw != "n":
                    repeattimes = repeattimes_raw.split(":")
                    try:
                        repeattimes[0] = int(repeattimes[0])
                    except:
                        repeattimes[0] = None

                    try:
                        repeattimes[1] = int(repeattimes[1])
                    except:
                        repeattimes[1] = None

                repeat_min = repeattimes[0]
                repeat_max = repeattimes[1]

                if repeat_min == None:
                    repeat_min = 1
                if repeat_max == None:
                    repeat_max = 2e100

                if repeat_min > repeat_max:
                    raise BadRepeatTimesError(li)

                parentmono = mono[end_ind]
                childmono = mono[start_ind]

                l = None
                if wurcssubststr == None:
                    l = parentmono.add_child_repeat_bridge(childmono,
                                         child_pos=start_pos,
                                         parent_pos=end_pos,
                                         parent_type=Linkage.oxygenPreserved,
                                         child_type=Linkage.oxygenLost)
                else:
                    wurcssubststr = wurcssubststr.replace("*", "")
                    subst = self.mf.getsubst(wurcssubststr)

                    substparenttype1 = eval(self.mf.subsconfig.get(wurcssubststr, "parent_type"))
                    substchildtype1 = eval(self.mf.subsconfig.get(wurcssubststr, "child_type"))

                    substparenttype2 = Linkage.nitrogenAdded
                    substchildtype2 = Linkage.oxygenPreserved

                    parentmono.add_substituent(subst,
                                               parent_type=substparenttype1,
                                               parent_pos=end_pos,
                                               child_type=substchildtype1,
                                               child_pos=1
                                               )
                    l = subst.add_child_repeat_bridge(childmono,
                                                      parent_type=substparenttype2,
                                                      parent_pos=1,
                                                      child_type=substchildtype2,
                                                      child_pos=start_pos)

                l.set_repeat_times_min(repeattimes[0])
                l.set_repeat_times_max(repeattimes[1])

                rep.append((childmono, parentmono, start_pos, end_pos, repeattimes))

                continue

            mi = self.substbridgelinkre.search(li)
            if mi:
                #print mi.group()
                ind1 = self.char2int[mi.group(1)]
                pos1 = (int(mi.group(2)) if mi.group(2) != "?" else None)
                ind2 = self.char2int[mi.group(3)]
                pos2 = (int(mi.group(4)) if mi.group(4) != "?" else None)

                wurcssubststr = mi.group(5).replace("*", "")
                subst = self.mf.getsubst(wurcssubststr)

                substparenttype1 = eval(self.mf.subsconfig.get(wurcssubststr, "parent_type"))
                substchildtype1 = eval(self.mf.subsconfig.get(wurcssubststr, "child_type"))
                substparenttype2 = Linkage.nitrogenAdded
                substchildtype2 = Linkage.oxygenPreserved

                if not (ind1 < ind2):
                    raise MonoOrderLinkError(li)

                parentmono = mono[ind1]
                childmono = mono[ind2]

                parentmono.add_substituent(subst, parent_type=substparenttype1, parent_pos=pos1, child_type=substchildtype1, child_pos=1)
                subst.add_child(childmono, parent_type=substparenttype2, parent_pos=1, child_type=substchildtype2, child_pos=pos2)

                continue

            raise UnsupportedLinkError(li)

        if counts[2] != (counts[1] - 1 + len(floating_substs) + len(rep)):
            raise LinkCountError(s)

        unconnected = set()
        monocnt = 0
        for id,m in mono.items():
            monocnt += 1
            # TODO not reliable
            if len(m.parent_links()) > 0 and not (len(m.parent_links()) == 1 and m.parent_links()[0].repeat_bridge_link()):
                continue
            m.set_connected(False)
            unconnected.add(m)

        if len(unconnected) not in (1,monocnt):
            tofix = set()
            for m in unconnected:
                if len(m.links()) != 1:
                    continue
                for l in m.links():
                    if l.parent_pos() == None and l.child_pos() == None:
                        tofix.add(m)
            
            if len(tofix) == (len(unconnected)-1):
                for m in tofix:
                    m.links()[0].reverse()
                    m.set_connected(True)
                    unconnected.remove(m)

        # print len(unconnected)
        # for m in unconnected:
        #     print m

        if len(unconnected) == 0:
            pass
            # TODO circular vs repeat differences?
            raise CircularError(s)
        # assert len(unconnected) in (1,monocnt), "# of unconnected nodes: %s not in {1,%d}"%(len(unconnected),monocnt)
        if len(unconnected) not in (0,1,monocnt):
            raise UnexpectedConnectivityError(s)
            
        if len(unconnected) == 1:
            g = Glycan(unconnected.pop())
            if len(floating_substs) > 0:
                raise UnexpectedFloatingSubstError(s)
            g.set_undetermined(undets)
        else:
            assert len(undets) == 0
            g = Glycan()
            g.set_undetermined(set(list(unconnected)+floating_substs))

        if g.is_repeat_in_repeat():
            raise UnexpectedRepeatInRepeatError()

        return g

    def add_child(self, parent, child, parent_pos, child_pos, repeat_info):
        exit_link = False

        for rinfo in repeat_info:

            if parent == rinfo[1] and parent_pos == rinfo[3] :
                # TODO check for child pos? child_pos == rinfo[2] bad example G45957UA
                # TODO perhaps also check the monosaccharide ?
                exit_link = True


        if exit_link:
            parent.add_child_exit_repeat(child,
                             parent_pos=parent_pos,
                             child_pos=child_pos,
                             parent_type=Linkage.oxygenPreserved,
                             child_type=Linkage.oxygenLost
                             )
        else:
            parent.add_child(child, parent_pos = parent_pos, child_pos = child_pos, parent_type = Linkage.oxygenPreserved, child_type = Linkage.oxygenLost)
        return

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
            print("+++", os.path.split(f)[1])
            # for t in g.undetermined_root_reprs():
            #     print t[1],str(t[0])
            print(GlycoCTFormat().toStr(g))
            for m in g.all_nodes(undet_subst=True):
                print(m)
            print(g.underivitized_molecular_weight())
        except GlycanParseError as e:
            print("!!!", os.path.split(f)[1], e)
            bad += 1
    print("Failed: %d/%d"%(bad,len(files)))
