
import sys, re

from .GlycanFormatter import *
from .GlycanBuilderSVGParser import *
from .CompositionFormatter import *
from .GlycanResource import GlyTouCanNoPrefetch as GlyTouCan

class GlycanMultiParser(object):
    wp = WURCS20Format()
    gp = GlycoCTFormat()
    cp = CompositionFormat()
    ip1 = IUPACLinearFormat()
    ip2 = IUPACParserExtended1()
    ip3 = IUPACParserGlyTouCanExtended()
    ip4 = IUPACParserCFG()
    svgp = GlycanBuilderSVG()
    gtc = GlyTouCan()

    def __init__(self):
        self.lastparser = None
    
    def parser(self):
        if self.lastparser is not None:
            return self.lastparser
        return None
    
    def parser_name(self):
        if self.lastparser is not None:
            return self.lastparser.__class__.__name__
        return None
    
    def toGlycan(self,seq):
        if seq.startswith("RES"):
            self.lastparser = self.gp
            return self.gp.toGlycan(seq)
        elif seq.startswith("WURCS"):
            self.lastparser = self.wp
            return self.wp.toGlycan(seq)
        elif seq.startswith("<svg"):
            self.lastparser = self.svgp
            return self.svgp.toGlycan(seq)
        elif re.search(r'^G\d{5}\w{2}$',seq):
            seq = self.gtc.getseq(seq,'wurcs')
            self.lastparser = self.wp
            return self.wp.toGlycan(seq)
        try:
            self.lastparser = self.cp
            return self.cp.toGlycan(seq)
        except GlycanParseError:
            pass
        try:
            self.lastparser = self.ip1
            return self.ip1.toGlycan(seq)
        except GlycanParseError:
            pass
        try:
            self.lastparser = self.ip4
            return self.ip4.toGlycan(seq)
        except GlycanParseError:
            pass
        try:
            self.lastparser = self.ip2
            return self.ip2.toGlycan(seq)
        except GlycanParseError:
            pass
        try:
            self.lastparser = self.ip3
            return self.ip3.toGlycan(seq)
        except GlycanParseError:
            pass
        self.lastparser = None
        raise GlycanParseError

    def normalizedSequence(self,seq):
        if seq.startswith("RES"):
            if '\\n' in seq:
                seq = seq.replace('\\n','\n')
            return "\n".join(seq.split())
        elif seq.startswith("WURCS"):
            return seq
        elif re.search(r'^G\d{5}\w{2}$',seq):
            return self.gtc.getseq(seq,'wurcs')
        try:
            return self.cp.toSequence(seq)
        except GlycanParseError:
            pass
        gly = self.toGlycan(seq)
        return gly.glycoct()
