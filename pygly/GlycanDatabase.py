
from GlycoCTDatabase import *
from CFGArrayDatabase import CFGArrayDatabase

__all__ = [ "GlycanDatabaseFactory" ];

extn2cls = {
             'gct': GlycoCTDatabase,
             'gdb': GlycomeDBDatabase,
	     'uca': UniCarbKBDatabase,
	     'cfg': CFGArrayDatabase,
	   }

def GlycanDatabaseFactory(filename):
    extn = filename.rsplit('.',1)[-1]
    if extn in extn2cls:
	return extn2cls[extn](filename)
    return None

