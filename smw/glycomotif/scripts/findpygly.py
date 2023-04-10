
try:
    import pygly
except ImportError:
    import os,os.path,sys
    base = os.path.split(os.path.abspath(sys.argv[0]))[0]
    while not os.path.exists(os.path.join(base,"pygly")):
        base = os.path.split(base)[0]
        if base in ("/",""):
            break
    if os.path.exists(os.path.join(base,"pygly")):
        sys.path.append(base)
    import pygly

