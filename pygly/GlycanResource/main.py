
from __future__ import print_function

import os.path

filename = os.path.join(os.path.dirname(__file__),"__main__.py")
with open(filename) as f:
    code = compile(f.read(), filename, 'exec')
    exec(code)
