#!/bin/env python27
import sys
from pygly.GlycoCTDatabase import GlycoCTDatabase, GlycomeDBDatabase
from pygly.CartoonistDatabase import CartoonistDatabase
from pygly.LinearIUPACDatabase import CFGArrayDatabase, GlycoSuiteDatabase
from pygly.GlycanFormatter import LinearCodeFormat, GlycoCTFormat
from pygly.GlyDbIndex import GlyDbIndex
from pygly.GlycosidaseFilter import NeuraminidaseFilter,GalactosidaseFilter,FucosidaseFilter,Beta14GalactosidaseFilter
from zipfile import ZipFile, ZIP_DEFLATED

glycosidase_options = filter(None,"""
neuraminidase       -Neuraminidase       NeuraminidaseFilter
sialidase           -Sialidase           NeuraminidaseFilter
galactosidase       -Galactosidase       GalactosidaseFilter
fucosidase          -Fucosidase          FucosidaseFilter
beta14galactosidase -Beta14Galactosidase Beta14GalactosidaseFilter
""".splitlines())

glyopts = {}
for opt in glycosidase_options:
    sopt = opt.split()
    glyopts[sopt[0]] = (sopt[1],eval(sopt[2]))
glycosidase_options =glyopts

if len(sys.argv) < 2:
    print >>sys.stderr, "No command, use one of: build, force-build, dump, count, glycoct."
    sys.exit(1)

if sys.argv[1] == "build" or sys.argv[1] == "force-build":
    force = (sys.argv[1] == "force-build")
    sys.argv.pop(1)
    filters = None
    modifier = None
    if len(sys.argv) > 1 and sys.argv[1] in glycosidase_options:
        filters = [glycosidase_options[sys.argv[1]][1]]
        modifier = glycosidase_options[sys.argv[1]][0]
        sys.argv.pop(1)
    if len(sys.argv) <= 1:
	print >>sys.stderr, "GlyDbIndex.py build [galactosidase-option] <database>.(gct|gdb|cfg)"
        print >>sys.stderr, "\nGlycosidase options:\n  "+'\n  '.join(sorted(glycosidase_options))
	sys.exit(1)
    if sys.argv[1].endswith('.gct'):
        gdb = GlycoCTDatabase(sys.argv[1])
    elif sys.argv[1].endswith('.gdb'):
        gdb = GlycomeDBDatabase(sys.argv[1])
    elif sys.argv[1].endswith('.cfg'):
        gdb = CFGArrayDatabase(sys.argv[1])
    elif sys.argv[1].endswith('.gsdb'):
        gdb = GlycoSuiteDatabase(sys.argv[1])
    else:
	if len(sys.argv) <= 2:
	    print >>sys.stderr, "GlyDbIndex.py build <database>.ctn <max-relative-demerits>"
	    sys.exit(1)
        gdb = CartoonistDatabase(sys.argv[1],maxreldem=int(sys.argv[2]))
    if modifier:
	sfn = gdb.filename.rsplit('.',1)
	indexfile = sfn[0] + modifier + '.' + sfn[1] + '.index'
        gdb = GlyDbIndex(gdb,force=force,filters=filters,indexfile=indexfile)
    else:
        gdb = GlyDbIndex(gdb,force=force,filters=filters)
    # gdb.build()
elif sys.argv[1] in ("dump","count","glycoct"):
    dump = False
    count = False
    glycoct = False
    if sys.argv[1] == "dump":
	dump = True
    if sys.argv[1] == "count":
	count = True
    if sys.argv[1] == "glycoct":
	glycoct = True
    if len(sys.argv) <= 1:
	print >>sys.stderr, "GlyDbIndex.py %s <glyco-database>"%(sys.argv[1],)
    sys.argv.pop(1)
    dbname = sys.argv[1]
    gdb = GlyDbIndex(sys.argv[1])
    sys.argv.pop(1)

    getargs = {}
    for i in range(1,len(sys.argv),2):
        key = sys.argv[i]
	try:
            value = sys.argv[i+1]
            value = float(sys.argv[i+1])
	    value = int(sys.argv[i+1])
	except ValueError:
	    pass
        getargs[key] = value
    if count:
        cnt = gdb.count(**getargs)
	print dbname,cnt
    elif glycoct:
        fmt1 = GlycoCTFormat()
	zf = ZipFile(dbname.rsplit('.',1)[0]+'.gct','w',ZIP_DEFLATED)
        for r in gdb.get(**getargs):
          zf.writestr("%s.txt"%r.accession,fmt1.toStr(r.glycan))
	zf.close()
    else:
        fmt = LinearCodeFormat()
        for r in gdb.get(**getargs):
          if r['lincode']:
            lc = fmt.toStr(r.glycan)
            print r.accession,r['nlinked'],r.get('oxford',"-"),r['molecular_weight'],r['composition'],lc
	    print r
          else:
            print r.accession,r['nlinked'],r.get('oxford',"-"),r['molecular_weight'],r['composition']
	    print r
else:
    print >>sys.stderr, "Bad command:"+sys.argv[1]
    sys.exit(1)
