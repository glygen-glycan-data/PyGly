
class GlycanParseError(Exception):
    def __str__(self):
        return self.message

class GlycoCTParseError(GlycanParseError):
    pass

class LinearCodeParseError(GlycanParseError):
    pass

class IUPACLinearParseError(GlycanParseError):
    pass

class IUPACParseError(GlycanParseError):
    pass

class WURCS20ParseError(GlycanParseError):
    pass

class GlycanBuilderSVGParseError(GlycanParseError):
    pass

