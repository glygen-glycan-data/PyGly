#!/bin/env python2

import sys,urllib,os,os.path,subprocess
import csv
import zipfile

zf = zipfile.ZipFile(sys.argv[1],'w')
header = None
for r in sys.stdin:
    r=list(map(lambda s: s.strip('"'),r.strip('\n').split('\t')))
    if not header:
        header = list(r)
        sys.stdout.write(header[0].lstrip('?')+"\n")
        continue
    zf.writestr("%s.txt"%(r[0],),eval('"' + r[1] + '"').strip()+"\n")
    sys.stdout.write(r[0]+"\n")
zf.close()
