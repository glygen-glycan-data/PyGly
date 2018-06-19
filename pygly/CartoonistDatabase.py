
import zipfile
import traceback
from operator import itemgetter

from GlycanFormatter import LinearCodeFormat, LinearCodeParseError
from GlyRecord import GlyRecord

class CartoonistDatabase:
    def __init__(self,filename,maxdem=9999,maxreldem=9999):
        self.filename = filename
        self.name,extn = filename.rsplit('.',1)
        assert extn == 'ctn'        
        self.fmt = LinearCodeFormat()
	self.maxdem = maxdem
	self.maxreldem = maxreldem
    def __iter__(self):
        return self.next()
    def group(self):
	atcomp = []
	curcomp = None
        h = open(self.filename)
        for lineno,l in enumerate(h):
	    if lineno == 0:
		continue
	    sl = l.split()
	    lc = sl[0]; mass = float(sl[1]); dem = int(sl[2])
	    comp = tuple([mass]+map(int,sl[3:]))
	    if comp != curcomp:
		if curcomp != None:
		    for r in self.process(atcomp):
			yield r
		curcomp = comp
		atcomp = [(lineno+1,lc,dem)]
	    else:
		atcomp.append((lineno+1,lc,dem))
	if len(atcomp) > 0:
	    for r in self.process(atcomp):
		yield r
    def process(self,rows):
	rows.sort(key=itemgetter(2))
	mindem = rows[0][2]
	for r in rows:
	    yield (r[0],r[1],r[2],r[2]-mindem)
    def next(self):
	for r in self.group():
	    acc,lc,dem,reldem = r
	    if dem > self.maxdem:
		continue
	    if reldem > self.maxreldem:
		continue
            try:
                g = self.fmt.toGlycan(lc)
	    except LinearCodeParseError:
                # traceback.print_exc()
		continue
            yield GlyRecord(source='Cartoonist',accession="CTN%d"%acc,glycan=g,demerits=dem,relative_demerits=reldem,name=self.name)

if __name__ == '__main__':
    import sys
    from GlycanFormatter import LinearCodeFormat
    fmt = LinearCodeFormat()
    gdb = CartoonistDatabase(sys.argv[1],maxreldem=0,maxdem=9999)
    for r in gdb:
        lc = fmt.toStr(r.glycan)
        # print r.accession,r['demerits'],r['relative_demerits'],lc
	# print r.glycan
	sys.stdout.flush()
