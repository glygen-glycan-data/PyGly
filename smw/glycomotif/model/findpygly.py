
import os,os.path,sys
pyglybase = os.path.split(os.path.abspath(sys.argv[0]))[0]
while not os.path.exists(os.path.join(pyglybase,"pygly")):
    pyglybase = os.path.split(pyglybase)[0]
    if pyglybase in ("/",""):
        break
if os.path.exists(os.path.join(pyglybase,"pygly")):
    sys.path.append(pyglybase)
