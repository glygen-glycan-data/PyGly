
from pkg_resources import resource_stream
from ConfigParser import SafeConfigParser
import os.path, sys

class ReferenceTable(dict):
    def __init__(self,iniFile=None):
        if not iniFile:
            iniFile = resource_stream(__name__, self.__class__.__name__.lower()+'.ini')
        elif isinstance(iniFile, basestring) and os.path.exists(iniFile):
            iniFile = open(iniFile)
        cfg = SafeConfigParser()
        cfg.optionxform = str
        cfg.readfp(iniFile)
        self.parseConfig(cfg)
    def parseConfig(self,cfg):
        for name in cfg.sections():
            self.update(dict(self.parseSection(name,dict(cfg.items(name)))))
