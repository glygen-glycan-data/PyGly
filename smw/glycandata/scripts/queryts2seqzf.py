#!/bin/env python2

import sys,urllib,os,os.path,subprocess
import csv
import zipfile

reader = csv.DictReader(sys.stdin,dialect='excel-tab')
zf = zipfile.ZipFile(sys.argv[1],'w')
header = None
for r in reader:
    if not header:
        header = list(reader.fieldnames)
        sys.stdout.write(header[0].lstrip('?')+"\n")
    if r.get(header[1]) != "":
	zf.writestr("%s.txt"%(r[header[0]],),eval('"'+r[header[1]]+'"').strip()+"\n")
    sys.stdout.write(r[header[0]]+"\n")
zf.close()
