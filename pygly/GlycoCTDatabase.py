
import zipfile, csv, copy
from cStringIO import StringIO
import traceback
import sys, os, os.path

from GlycanFormatter import GlycoCTFormat, GlycoCTParseError
from GlyRecord import GlyRecord
from CompositionTable import ResidueCompositionTable
from ElementMass import MonoisotopicElementMass
from GlyMWFilter import GlyMWFilter
from GlyCompFilter import GlyCompFilter
from GlyNLinkedFilter import GlyNLinkedFilter

# import time###
# f = open('parsing_errors.txt','w')

class GlycoCTDatabase:
    prefix = ""
    extn = "gct"
    source = ""
    def __init__(self,filename):
        self.filename = filename
        self.name,extn = filename.rsplit('.',1)
        self.name = os.path.split(self.name)[1]
        assert extn == self.extn
        self.fmt = GlycoCTFormat()
    def getraw(self,accession):
	zf = zipfile.ZipFile(self.filename, "r")
	fn = accession[len(self.prefix):]+'.txt'
	try:
            zf.getinfo(fn)
        except KeyError:
            return None
	gct = self._getraw(zf,fn)
	zf.close()
	return gct
    def _getraw(self,zf,name):
         return zf.read(name)
    def _get(self,zf,name):
        try:
            glystr = self._getraw(zf,name)
            g = self.fmt.toGlycan(glystr)
        except GlycoCTParseError:
            # f.write('> GlycanFormatter > GLYCANPARSEERROR > '+self.message)###
            # print '> GlycanFormatter > GLYCANPARSEERROR > ',self.message###
            # time.sleep(2)###
            return None
        except:
            print >>sys.stderr, "Problem with GlycoCT file "+name
            # time.sleep(2)###
            traceback.print_exc()
            sys.exit(1)
	kwargs = {}
	try:
	    attstr = zf.read(name.rsplit('.',1)[0]+'.att')
	    for r in csv.reader(StringIO(attstr)):
	        kwargs[r[0]] = copy.copy(r[1:])
	except KeyError:
	    pass
	try:
	    kwargs['image'] = zf.read(name.rsplit('.',1)[0]+'.png')
	except KeyError:
	    pass
        gr = GlyRecord(source=self.source if self.source else self.name,
                       accession=self.prefix+name.rsplit('.',1)[0],glycan=g,name=self.name,
		       **kwargs)
	# print gr
	return gr
    def get(self,accession):
	zf = zipfile.ZipFile(self.filename, "r")
	fn = accession[len(self.prefix):]+'.txt'
	try:
	    zf.getinfo(fn)
	except KeyError:
	    return None
	gr = self._get(zf,fn)
	zf.close()
	return gr
    def __iter__(self):
        return self.next()
    def next(self):
        zf = zipfile.ZipFile(self.filename, "r")
        for name in zf.namelist():
            if not name.endswith('.txt'):
                continue
	    gr = self._get(zf,name)
	    if gr:
	 	yield gr
	zf.close()

class GlycomeDBDatabase(GlycoCTDatabase):
    prefix = "GDB"
    extn = "gdb"
    source = "GlycomeDB"

class UniCarbKBDatabase(GlycoCTDatabase):
    prefix = "UC"
    extn = "uca"
    source = "UniCarbKB"

if __name__ == '__main__':
    import sys
    from GlycanFormatter import LinearCodeFormat
    compkeys = ('HexNAc','Hex','Fuc','NeuAc','NeuGc','Xyl','Xxx')
    fmt = LinearCodeFormat()
    extn = sys.argv[1].rsplit('.',1)[-1]
    if extn == 'gdb':
        gdb = GlycomeDBDatabase(sys.argv[1])
    elif extn == 'gct':
        gdb = GlycoCTDatabase(sys.argv[1])
    else:
	raise RuntimeError('Unknown extension...')
    gdb = GlyMWFilter(gdb,ResidueCompositionTable(),MonoisotopicElementMass(),addh2o=True)
    gdb = GlyCompFilter(gdb)
    gdb = GlyNLinkedFilter(gdb)
    for r in gdb:
        # if not r['nlinked']:
        #     continue
        try:
            lc = fmt.toStr(r.glycan)
            print r.accession,r['nlinked'],r['molecular_weight'],r['composition'].str(compkeys),lc
	    # print r.glycan
        except KeyError:
            print r.accession,r['nlinked'],r['molecular_weight'],r['composition'].str(compkeys),None
