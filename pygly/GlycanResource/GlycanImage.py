
try:
    from pygly.GlycanImage import GlycanImage as GI
except:
    # from GlycanImage import GlycanImage as GI
    pass

import tempfile, os.path, traceback
from .GlyTouCan import GlyTouCan

class GlycanImage(object):
    def __init__(self):
        self.imageWriter = GI()
        self.gtc = GlyTouCan(prefetch=False)
    def write(self,*args,**kw):
        kw['format'] = kw.get('format','png')
        kw['force'] = kw.get('force',True)
	out = None
	if 'out' in kw:
	    out = kw.get('out')
            del kw['out']
        for k,v in list(kw.items()):
            self.imageWriter.set(k,v)
        written = []
        for acc in args:
            if not out:
               outfile = "%s.%s"%(acc,kw['format'])
	    else:
               outfile = out

            seq = self.gtc.getseq(acc,'wurcs')
            if seq:
                self.imageWriter.writeImage(seq,outfile)
                if os.path.exists(outfile) and os.path.getsize(outfile) > 0:
                    written.append(outfile)
                    continue

            seq = self.gtc.getseq(acc,'glycoct')
            if seq:
                self.imageWriter.writeImage(seq,outfile)
                if os.path.exists(outfile) and os.path.getsize(outfile) > 0:
                    written.append(outfile)
                    continue

            seq = self.gtc.glycoct(acc)
            if seq:
                self.imageWriter.writeImage(seq,outfile)
                if os.path.exists(outfile) and os.path.getsize(outfile) > 0:
                    written.append(outfile)
                    continue

        return written

    def show(self,accession,**kw):

        kw['format'] = kw.get('format','png')
        dummy,tmpoutfile = tempfile.mkstemp(suffix=".%s"%(kw['format'],))
        kw['out'] = tmpoutfile
        kw['force'] = True
	args = [ accession ]
        self.write(*args,**kw)

        if not os.path.exists(tmpoutfile) or os.path.getsize(tmpoutfile) == 0:
            os.unlink(tmpoutfile)
	    return ""

        def button_click_exit_mainloop (event):
            event.widget.quit() # this will cause mainloop to unblock.

        import Tkinter
        from PIL import Image, ImageTk

        try:
            root = Tkinter.Tk()
            root.bind("<Button>", button_click_exit_mainloop)
            root.geometry('+%d+%d' % (100,100))
            image1 = Image.open(tmpoutfile)
            root.geometry('%dx%d' % (image1.size[0],image1.size[1]))
            tkpi = ImageTk.PhotoImage(image1)
            label_image = Tkinter.Label(root, image=tkpi)
            label_image.place(x=0,y=0,width=image1.size[0],height=image1.size[1])
            root.title(accession)
            root.mainloop() # wait until user clicks the window
	finally:
            os.unlink(tmpoutfile)

	return ""
