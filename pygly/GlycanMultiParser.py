
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
    svgp = GlycanBuilderSVG()
    gtc = GlyTouCan()

    def toGlycan(self,seq):
        if seq.startswith("RES"):
            return self.gp.toGlycan(seq)
        elif seq.startswith("WURCS"):
            return self.wp.toGlycan(seq)
        elif seq.startswith("<svg"):
            return self.svgp.toGlycan(seq)
        elif re.search(r'^G\d{5}\w{2}$',seq):
            seq = self.gtc.getseq(seq,'wurcs')
            return self.wp.toGlycan(seq)
        try:
            return self.cp.toGlycan(seq)
        except GlycanParseError:
            pass
        try:
            return self.ip1.toGlycan(seq)
        except GlycanParseError:
            pass
        try:
            return self.ip2.toGlycan(seq)
        except GlycanParseError:
            pass
        raise GlycanParseError

    def normalizedSequence(self,seq):
        if seq.startswith("RES"):
            return seq
        elif seq.startswith("WURCS"):
            return seq
        try:
            return self.cp.toSequence(seq)
        except GlycanParseError:
            pass
        gly = self.toGlycan(seq)
        return gly.glycoct()
