
from __future__ import print_function

from . GlycanResource import GlyTouCanNoCache
from . GlycanFormatter import GlycanFormatter, WURCS20Format
from . GlycanFormatterExceptions import *
from . Glycan import Glycan

import re, sys

class BadComposition(CompositionParseError):
    def __init__(self,instr):
        self.message = "Composition parser: Bad composition: %s"%(instr,)

class CompositionFormat(GlycanFormatter):
    def __init__(self):
        self.wp = WURCS20Format()
        self.gtc = GlyTouCanNoCache()

    def toStr(self,composition=None,glycan=None,**kwargs):
        if not composition:
            composition = glycan.iupac_composition(**kwargs)
        keys = set([k for k in composition.keys() if composition[k] > 0])
        keys = keys-set(["Count"])
        complist = []
        seen = set()
        for sym in Glycan.iupac_composition_syms + \
                   Glycan.iupac_aldi_composition_syms + \
                   Glycan.subst_composition_syms + \
                   [ "Xxx", "X" ]:
            if composition[sym] > 0:
                complist.append("%s(%d)"%(sym,composition[sym]))
                seen.add(sym)
        # print(seen,keys,complist)
        assert seen == keys
        return "".join(complist)
        
    def toSequence(self,s):
        if s.startswith('comp_'):
            s = s[5:]
        if '(' in s:
            sl = re.split(r'\s*\((\s*\d+\s*)\)\s*',s)
        else:
            sl = re.split(r'\s*(\d+)\s*',s)
        sl = list(map(str.strip,sl))
        if len(sl)%2 != 1 or sl[-1].strip() != "":
            raise BadComposition(s)
        comp = dict()
        for i in range(0,len(sl)-1,2):
            try:
                cnt = int(sl[i+1])
            except ValueError:
                raise BadComposition(s)
            if cnt < 0:
                raise BadComposition(s)
            if cnt == 0:
                continue
            comp[sl[i]] = cnt
        try:
            seq = self.gtc.makecompseq(**comp)
        except (ValueError, TypeError, KeyError):
            raise BadComposition(s)
        return seq

    def toGlycan(self,s):
        seq = self.toSequence(s)
        try:
            gly = self.wp.toGlycan(seq)
        except GlycanParseError:
            raise BadComposition(s)
        return gly

class FiveLetterComposition(CompositionFormat):
    regex = re.compile(r'^N(\d+)H(\d+)F(\d+)S(\d+)G(\d+)$')
    
    def toSequence(self,s):
        if not self.regex.search(s.strip()):
            raise BadComposition(s)
        m = self.regex.search(s.strip())
        s1 = f"HexNAc{m.group(1)}Hex{m.group(2)}Fuc{m.group(3)}NeuAc{m.group(4)}NeuGc{m.group(5)}"
        return super().toSequence(s1)
    
    def toStr(self,*args,**kwargs):
        compstr = super().toStr(*args,**kwargs)
        splcomp = re.split(r'\((\d+)\)',compstr)
        comp = dict()
        for i in range(0,len(splcomp)-1,2):
            cnt = int(sl[i+1])
            if cnt == 0:
                continue
            comp[sl[i]] = str(cnt)
        if set(comp.keys()) <= set(['HexNAc','Hex','Fuc','NeuAc','NeuGc']):
            return f"N{comp['HexNAc']}H{comp['Hex']}F{comp['Fuc']}S{comp['NeuAc']}G{comp['NeuGc']}"
        raise BadComposition(compstr)