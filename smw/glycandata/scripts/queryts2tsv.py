#!/bin/env python3.12

import sys
import csv

reader = csv.DictReader(sys.stdin,dialect='excel-tab')
header = None
for r in reader:
    if not header:
        header = map(lambda h: h.lstrip('?'),reader.fieldnames)
        sys.stdout.write("\t".join(header)+"\n") 
    line = map(lambda h: r[h].strip('"'),reader.fieldnames)
    sys.stdout.write("\t".join(line)+"\n")
