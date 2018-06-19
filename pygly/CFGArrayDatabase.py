
import zipfile
import traceback
from operator import itemgetter
import re

from GlycanFormatter import IUPACLinearFormat, IUPACLinearParseError, GlycoCTFormat
from GlyRecord import GlyRecord

class CFGArrayDatabase:
    def __init__(self,filename):
        self.filename = filename
        self.name,extn = filename.rsplit('.',1)
        assert extn == 'cfg'
        self.fmt = IUPACLinearFormat()
    def __iter__(self):
        return self.next()
    def next(self):
        h = open(self.filename)
        for lineno,l in enumerate(h):
	    if lineno == 0:
		continue
	    if l.startswith('#'):
		continue
	    sl = l.split(None,1)
	    acc = int(sl[0]); glystr = sl[1].strip()
            # Wow - badly hand formatted CFG IUPAC linear glycans suck...
            m = re.search(r'([ab]1?)?-?(\d(-?))?([RS][pP]?\d+|MDPLys)?$',glystr)
            if m != None:
                x = len(m.group(0))
                glystr = glystr[:-x]
            try:
                g = self.fmt.toGlycan(glystr)
	    except IUPACLinearParseError:
                traceback.print_exc()
                continue
            print >>sys.stderr, ">>CFG%03d"%acc
            yield GlyRecord(source="CFGArray",accession="CFG%03d"%acc,glycan=g,name=self.name)            

if __name__ == '__main__':
    import sys
    from GlycanFormatter import IUPACLinearFormat
    fmt = IUPACLinearFormat()
    fmt1 = GlycoCTFormat()
    gdb = CFGArrayDatabase(sys.argv[1])
    for r in gdb:
        lc = fmt.toStr(r.glycan)
        print r.accession,lc
	# print fmt1.toStr(r.glycan)
	# print r.glycan
	print r.glycan
	sys.stdout.flush()
