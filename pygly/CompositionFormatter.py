
from __future__ import print_function

from . GlycanResource import GlyTouCanNoCache
from . GlycanFormatter import GlycanFormatter, WURCS20Format
from . GlycanFormatterExceptions import *

import re, sys

class BadComposition(CompositionParseError):
    def __init__(self,instr):
        self.message = "Composition parser: Bad composition: %s"%(instr,)

class CompositionFormat(GlycanFormatter):
    def __init__(self):
        self.wp = WURCS20Format()
        self.gtc = GlyTouCanNoCache()
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

