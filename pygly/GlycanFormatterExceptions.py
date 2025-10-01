
class GlycanParseError(Exception):
    def __init__(self,*args,**kwargs):
        self.message = "Glycan parse error"
        super().__init__(*args,**kwargs)

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

class CompositionParseError(GlycanParseError):
    pass
