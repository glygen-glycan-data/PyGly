
import re
from collections import defaultdict
from ConstantsTable import ConstantsTable

class Composition(defaultdict):
    def __init__(self,*args,**kw):
        super(Composition,self).__init__(*args,**kw)
        self.default_factory = int
    def parse(self,line):
        sl= line.split()
        for i in range(0,len(sl),2):
            count = int(sl[i+1])
            if count == 0:
                continue
            self[sl[i]] += count
    def compactparse(self,word):
	sl = map(str.strip,re.split(r'([A-Z][a-z]?)',word)[1:])
	for i in range(1,len(sl),2):
	    if not sl[i]:
		sl[i] = '1'
	self.parse((" ".join(sl)).strip())
    @staticmethod
    def fromstr(*args):
	c = Composition()
	for w in args:
	    c.compactparse(w)
	return c 
    def __str__(self):
        return ' '.join("%s %d"%(e,c) for e,c in sorted(self.items()) if c != 0)
    def str(self,keys):
        return ' '.join("%s %d"%(e,self[e]) for e in keys)
    def compactstr(self):
	return ''.join("%s%s"%(e,c if c > 1 else "") for e,c in sorted(self.items()) if c != 0)
    def add(self,c):
        for k in c:
            self[k] += c[k]
    def scale(self,s):
        for k in c:
            self[k] = int(s*self[k])
    def sub(self,c):
        for k in c:
            self[k] -= c[k]
    def contains(self,c):
        for k in c:
            if self[k] < c[k]:
                return False
        return True
    def eq(self,c):
        for k in c:
            if self[k] != c[k]:
                return False
        return True
    def mass(self,mass_table):
        return sum(mass_table[e]*c for e,c in self.iteritems())
    def count(self):
        return sum(c for e,c in self.iteritems())

class ResidueCompositionTable(dict):
    def __init__(self):
        consts = ConstantsTable()
        for sym,kv in consts.iteritems():
            if 'ResidueComposition' in kv:
                c = Composition()
                c.parse(kv['ResidueComposition'])
                self[sym] = c
    def new(self):
	return Composition()

class PermethylCompositionTable(dict):
    def __init__(self):
        consts = ConstantsTable()
        for sym,kv in consts.iteritems():
            if 'PermethylComposition' in kv and 'ResidueComposition' in kv:
                c = Composition()
                c.parse(kv['ResidueComposition'])
                c1 = Composition()
                c1.parse(kv['PermethylComposition'])
                c.add(c1)
                self[sym] = c
    def new(self):
	return Composition()
