#!/bin/env python27
import sys, os, tempfile, atexit, Tkinter
from PIL import Image, ImageTk
from pygly.GlyDbIndex import GlyDbIndex
from pygly.GlycanImage import GlycanImage
from pygly.GlycoCTDatabase import GlycoCTDatabase, GlycomeDBDatabase

tmpoutfile = None
def cleanup():
    if tmpoutfile:
        try:
            os.unlink(tmpoutfile)
        except:
	    pass

# atexit.register(cleanup)

if len(sys.argv) <= 2:
    print >>sys.stderr, "GlyDbImage.py [ --noindex ] <glyco-database> <accession> [ image options ]"
    print >>sys.stderr, "GlyDbImage.py --file <glycoct-file> [ image options ]"
    print >>sys.stderr, """
    Image options:
    format       (png|jpg|...)                         [png]
    scale        <float>                               [4.0]
    reducing_end (true|false)                          [true]
    orientation  (RL|LR|TB|BT)                         [RL]
    notation     (cfg|cfgbw|cfglink|uoxf|text|uoxfcol) [cfg]
    display      (normal|normalinfo|compact)           [normalinfo]
    out          <filename>                            [to screen]
    """.strip()
    sys.exit(1)

glystr = None
g = None
if sys.argv[1] == '--file':

    sys.argv.pop(1)
    acc = sys.argv[1]
    if not os.path.exists(acc):
	print >>sys.stderr, "Can't open file %s"%acc
        sys.exit(1)
    h = open(acc)
    glystr = h.read()
    h.close()
    sys.argv.pop(1)

elif sys.argv[1] == '--noindex':

    sys.argv.pop(1)
    if sys.argv[1].endswith('.gct'):
	gdb = GlycoCTDatabase(sys.argv[1])
    elif sys.sargv[1].endsiwth('.gdb'):
	gdb = GlycomeDBDatase(sys.argv[1])
    else:
	raise RuntimeError("Unsupported database type")
    sys.argv.pop(1)

    acc = sys.argv[1]
    sys.argv.pop(1)

    glystr = gdb.getraw(acc)
    if glystr == None:
        print >>sys.stderr, "Can't find glycan with accession %s"%acc
        sys.exit(1)

else:

    gdb = GlyDbIndex(sys.argv[1])
    sys.argv.pop(1)

    acc = sys.argv[1]
    sys.argv.pop(1)
    
    try:
        g = iter(gdb.get(accession=acc)).next()
    except StopIteration:
        print >>sys.stderr, "Can't find glycan with accession %s"%acc
        sys.exit(1)

outfile = None
imageWriter = GlycanImage()
for i in range(1,len(sys.argv),2):
    key = sys.argv[i]
    value = sys.argv[i+1]
    if key == "out":
	outfile = value
	continue
    imageWriter.set(key,value)
fmt = imageWriter.format()
if not outfile:
    dummy,tmpoutfile = tempfile.mkstemp(suffix=".%s"%fmt)
    outfile = tmpoutfile
imageWriter.writeImage(g.glycan if g else glystr,outfile)
if not tmpoutfile:
    sys.exit(1)

def button_click_exit_mainloop (event):
    event.widget.quit() # this will cause mainloop to unblock.

root = Tkinter.Tk()
root.bind("<Button>", button_click_exit_mainloop)
root.geometry('+%d+%d' % (100,100))
image1 = Image.open(tmpoutfile)
root.geometry('%dx%d' % (image1.size[0],image1.size[1]))
tkpi = ImageTk.PhotoImage(image1)
label_image = Tkinter.Label(root, image=tkpi)
label_image.place(x=0,y=0,width=image1.size[0],height=image1.size[1])
root.title(acc)
root.mainloop() # wait until user clicks the window
