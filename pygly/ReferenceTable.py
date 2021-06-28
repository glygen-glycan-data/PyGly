
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)
from pkg_resources import resource_stream
try:
    from configparser import ConfigParser, SafeConfigParser
except ImportError:
    from ConfigParser import ConfigParser, SafeConfigParser

try:
    from past.builtins import basestring
except ImportError:
    pass

import os.path, sys

try:
    from StringIO import StringIO
except ImportError:
    pass

class ReferenceTable(dict):
    def __init__(self,iniFile=None):
        if not iniFile:
            iniFile = [ s.decode('utf8') for s in resource_stream(__name__, self.__class__.__name__.lower()+'.ini').read().splitlines() ]
            self.iniFile = self.__class__.__name__.lower()+'.ini'
        elif isinstance(iniFile, basestring) and os.path.exists(iniFile):
            iniFile = open(iniFile)
            self.iniFile = iniFile
        cfg = SafeConfigParser()
        cfg.optionxform = str
        if hasattr(cfg,'read_file'):
            cfg.read_file(iniFile,self.iniFile)
        else:
            cfg.readfp(StringIO(u'\n'.join(iniFile)),self.iniFile)
        self.parseConfig(cfg)
    def parseConfig(self,cfg):
        for name in cfg.sections():
            for i,(k,v) in enumerate(self.parseSection(name,dict(cfg.items(name)))):
                type = ("section" if i == 0 else "alias")
                assert k not in self, "Repeated %s \"%s\" in %s"%(type,k,self.iniFile)
                self[k] = v
