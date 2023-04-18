#tiny script, cairosvg can be called directly from the command line
#but it requires python 3 and I didn't want to figure out how to deal with that 
#since we have multiple pythons

import sys
from cairosvg import svg2png

svg_file = sys.argv[1]
out_file = sys.argv[2]

svg2png(file_obj=open(svg_file, "rb"), write_to = out_file)


