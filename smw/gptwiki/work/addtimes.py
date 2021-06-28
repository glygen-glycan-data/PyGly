#!/bin/env python3
import sys
times = dict()
for l in sys.stdin:
    sl = l.rsplit(None,2)
    if sl[2] == "times":
        times[sl[0]] = times.get(sl[0],0)+int(sl[1])
    else:
        sys.stdout.write(l)
for k,v in sorted(times.items(),key=lambda t: (t[1],t[0])):
    print(k,v,"times")
