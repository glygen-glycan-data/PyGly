#!/bin/env python27

from getwiki import GlycoMotifWiki
import sys, os, os.path, zipfile

class GlyTouCanCache(object):
    dir =  ".gtccache"
    def __init__(self):
	self._id2gtc = {}
	if os.path.exists(os.path.join(self.dir,"glytoucan.tsv")):
	    for l in open(os.path.join(self.dir,"glytoucan.tsv")):
	    	id,gtc = l.split()
	        self._id2gtc[id] = gtc
	self._wurcszf = None
	if os.path.exists(os.path.join(self.dir,"wurcs.zip")):
	    self._wurcszf = zipfile.ZipFile(os.path.join(self.dir,"wurcs.zip"))
	self._glycoctzf = None
	if os.path.exists(os.path.join(self.dir,"glycoct.zip")):
	    self._glycoctzf = zipfile.ZipFile(os.path.join(self.dir,"glycoct.zip"))

    def id2gtc(self,id):
	return self._id2gtc.get(id)

    def gtc2wurcs(self,gtc):
        if not self._wurcszf:
	    return None
	try:
	    zi = self._wurcszf.getinfo("%s.txt"%(gtc,))
	except KeyError:
	    return None
	return self._wurcszf.read(zi)

    def gtc2glycoct(self,gtc):
        if not self._glycoctzf:
	    return None
	try:
	    zi = self._glycoctzf.getinfo("%s.txt"%(gtc,))
	except KeyError:
	    return None
	return self._glycoctzf.read(zi)

if __name__ == "__main__":

    # we must instantiate early, before any use of the command-line. 
    w = GlycoMotifWiki()
    #print w._prefix
    
    dir = ".gtccache"
    try:
        os.makedirs(dir)
    except OSError:
        pass
    
    wh = open(os.path.join(dir,"glytoucan.tsv"),'w')
    zf1 = zipfile.ZipFile(os.path.join(dir,"wurcs.zip"),'w')
    zf2 = zipfile.ZipFile(os.path.join(dir,"glycoct.zip"),'w')
    
    wurcsseen = set()
    glycoctseen = set()
    for m in w.itermotif():
        print >>sys.stderr, m.get('id')
        glytoucan = m.get('glytoucan')
        wh.write("%s\t%s\n"%(m.get('id'),glytoucan))
        if glytoucan not in wurcsseen:
            wurcs = m.get('wurcs')
            if wurcs:
	        zf1.writestr("%s.txt"%(glytoucan,),wurcs.strip())
	        wurcsseen.add(glytoucan)
        if glytoucan not in glycoctseen:
	    glycoct = m.get('glycoct')
	    if glycoct:
	        zf2.writestr("%s.txt"%(glytoucan,),glycoct.strip()+"\n")
	        glycoctseen.add(glycoct)
    
    zf1.close()
    zf2.close()
    wh.close()
