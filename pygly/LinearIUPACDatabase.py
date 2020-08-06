
import zipfile
import traceback
from operator import itemgetter
import re

from GlycanFormatter import IUPACLinearFormat, IUPACLinearParseError
from GlyRecord import GlyRecord

class LinearIUPACDatabase:
    extn = 'iupac'
    source = ""
    prefix = ""
    def __init__(self,filename):
        self.filename = filename
        self.name,extn = filename.rsplit('.',1)
        assert extn == self.extn
        self.fmt = IUPACLinearFormat()
    def __iter__(self):
        return self.next()
    def preprocess(self,glystr):
	return glystr
    def next(self):
        h = open(self.filename)
	seenaccs = set()
        for lineno,l in enumerate(h):
	    if lineno == 0:
		continue
	    if l.startswith('#'):
		continue
	    sl = [ s.strip() for s in l.split(None,1) ]
	    acc = self.acc(sl[0]); glystr = self.preprocess(sl[1])
	    if acc in seenaccs:
		continue
	    seenaccs.add(acc)
            try:
                g = self.fmt.toGlycan(glystr)
	    except IUPACLinearParseError:
                traceback.print_exc()
                continue
            yield GlyRecord(source=self.source,accession=self.prefix+acc,glycan=g,name=self.name)            
 	h.close()

class CFGArrayDatabase(LinearIUPACDatabase):
    extn = 'cfg'
    source = 'CFGArray'
    prefix = 'CFG'
    def preprocess(self,glystr):
        # Wow - badly hand formatted CFG IUPAC linear glycans suck...
        m = re.search(r'([ab]1?)?-?(\d(-?))?([RS][pP]?\d+|MDPLys)?$',glystr)
        if m != None:
            x = len(m.group(0))
            glystr = glystr[:-x]
	return glystr
    def acc(self,s):
	return "%03d"%int(s)

class IUPACArrayDatabase(LinearIUPACDatabase):
    extn = 'iupac'
    source = 'IUPACArray'
    prefix = ''
    def preprocess(self,glystr):
        m = re.search(r'(R1)?$',glystr)
        if m != None:
            x = len(m.group(0))
            glystr = glystr[:-x]
	return glystr
    def acc(self,s):
	return s

class GlycoSuiteDatabase(LinearIUPACDatabase):
    extn = 'gsdb'
    source = 'GlycoSuiteDB'
    prefix = 'GSD'
    def preprocess(self,structure):
        structure = re.sub(r'\(([?ab][?\d]-[?\d])\)','\g<1>',structure)
        structure = structure.replace('[','(')
        structure = structure.replace(']',')')
	return structure
    def acc(self,s):
	return s

class EuroCarbFragmentDesc(LinearIUPACDatabase):
    extn = 'frag'
    source = 'EuroCarbDB'
    prefix = 'ECDB'
    def preprocess(self,structure):
        print structure
	structure = structure.replace(',p','')
	toklist = []
	for tok in structure.split('--'):
	    stok = tok.split('-')
	    tok = '-'.join(reversed(stok))
	    toklist.append(tok)
	structure = '-'.join(reversed(toklist))
	# structure = '-'.join(reversed(structure.split('--')))
	print structure
	structure = re.sub(r'-\?b1D-freeEnd$','',structure)
	structure = re.sub(r'\(-([1-9?])b1',r'b\1-1)',structure)
	structure = re.sub(r'\(-([1-9?])a1',r'a\1-1)',structure)
	structure = re.sub(r'\)-([1-9?])b1',r'b\1-1(',structure)
	structure = re.sub(r'\)-([1-9?])a1',r'a\1-1(',structure)
	structure = re.sub(r'-([1-9?])b1',r'b\1-1',structure)
	structure = re.sub(r'-([1-9?])a1',r'a\1-1',structure)
	print structure
	structure = re.sub(r'D-','',structure)
	structure = re.sub(r'/#ycleavage','',structure)
	print structure
	return structure
    def acc(self,s):
	return s

if __name__ == '__main__':
    import sys
    from GlycanFormatter import IUPACLinearFormat
    fmt = IUPACLinearFormat()
    gdb = EuroCarbFragmentDesc(sys.argv[1])
    for r in gdb:
        lc = fmt.toStr(r.glycan)
        print r.accession,lc
	print r.glycan
	sys.stdout.flush()
