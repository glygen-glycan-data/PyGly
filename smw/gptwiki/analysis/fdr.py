#!/bin/env python27

from collections import defaultdict

class FDRComputation(object):
    _groupingkey = "peptide_group_label"                                                                              
    _rankingkey = "xx_lda_prelim_score"
    _rankingdir = -1
    _ndecoys = 1

    def __init__(self, **kw):
	self.groupingkey = kw.get('groupingkey',self._groupingkey)
	self.rankingkey = kw.get('rankingkey',self._rankingkey)
	self.rankingdir = (1 if (kw.get('rankingdir',self._rankingdir)) > 0 else -1)
	self.ndecoys = kw.get('ndecoys',self._ndecoys)
        self.clear()

    def clear(self):
	self.counts = defaultdict(lambda : defaultdict(int))
	self.scores = None
	self.cumulative_counts = None
	self.fdr = None

    def compute_cumulative_counts(self):
        self.scores = sorted(self.counts.keys(),reverse=(self.rankingdir<0))
        self.cumulative_counts = defaultdict(lambda : defaultdict(int))
	for decoy in (True,False):
            self.cumulative_counts[self.scores[0]][decoy] = self.counts[self.scores[0]][decoy]
        for i in range(1,len(self.scores)):
            for decoy in (True,False):
                self.cumulative_counts[self.scores[i]][decoy] = \
                                                self.cumulative_counts[self.scores[i-1]][decoy] + \
                                                self.counts[self.scores[i]][decoy]

    def compute_fdr(self):
	if self.cumulative_counts == None:
	    self.compute_cumulative_counts()
        self.fdr = dict()
        lastfdri = 1.0
        for sc in reversed(self.scores):
            fdri = (float(self.cumulative_counts[sc][True])/self.ndecoys)/float(self.cumulative_counts[sc][False])
            if fdri < lastfdri:
                self.fdr[sc] = fdri
                lastfdri = fdri

    def score(self,qv):
	if self.fdr == None:
	    self.compute_fdr()
	for sc in sorted(self.fdr,reverse=(self.rankingdir>0)):
	    if self.fdr[sc] < qv:
		return sc
	return max(self.fdr)

    def targets(self,x):
	if self.fdr == None:
	    self.compute_fdr()
	for sc in sorted(self.fdr,reverse=(self.rankingdir<0)):
	    if self.rankingdir*x <= self.rankingdir*sc:
		return self.cumulative_counts[sc][False]
	return 0.0

    def qvalue(self,x):
	if self.fdr == None:
	    self.compute_fdr()
	for sc in sorted(self.fdr,reverse=(self.rankingdir<0)):
	    if self.rankingdir*x < self.rankingdir*sc:
		return self.fdr[sc]
	return 1.0

    def allqvalues(self):
	if self.fdr == None:
	    self.compute_fdr()
	for sc in sorted(self.fdr,reverse=(self.rankingdir<0)):
	    yield sc,self.fdr[sc],self.cumulative_counts[sc][False],self.cumulative_counts[sc][True]

class SeparateAnalysisFDR(FDRComputation):

    def add_ids(self,rows,decoy=False):
	last_group = None
	for row in sorted(rows,key=lambda r: (r[self.groupingkey],self.rankingdir*float(r[self.rankingkey]))):
	    group = row[self.groupingkey]
            if group == last_group:
                continue
	    last_group = group
            score = float(row[self.rankingkey])
            self.counts[score][decoy] += 1
	    last_group = group

    def add_target_ids(self,rows):
	self.add_ids(rows,False)

    def add_decoy_ids(self,rows):
	self.add_ids(rows,True)

class CombinedAnalysisFDR(FDRComputation):

    def add_ids(self,rows):
	last_group = None
	for row in sorted(rows,key=lambda r: (r[self.groupingkey],self.rankingdir*float(r[self.rankingkey]))):
	    group = row[self.groupingkey]
            if group == last_group:
                continue
	    last_group = group
            score = float(row[self.rankingkey])
	    decoy = (int(row['decoy']) > 0)
            self.counts[score][decoy] += 1
	    last_group = group

if __name__ == "__main__":
    
    import sys
    import csv

    if len(sys.argv) < 2:
	print >>sys.stderr, """
Usage: fdr.py <target-ids> <decoy-ids-1> [ <decoy-ids-2> ... ]
       fdr.py <combined-ids> [ <ndecoys> ]

Use the first form if you have separate target and decoy analysis files,
and the second form if you have a single analysis file with both target and
decoy identifications.
""".strip()
        print >>sys.stderr,""
	sys.exit(1)

    # Assume the first file is target ids, and the remainder are decoy ids
    target_file = sys.argv[1]
    decoy_files = None
    try:
	ndecoys = int(sys.argv[2])
    except IndexError:
	ndecoys = 1
    except ValueError:
        decoy_files = sys.argv[2:]
        ndecoys = len(decoy_files)

    if decoy_files == None:
	# fdr = CombinedAnalysisFDR(ndecoys=ndecoys,rankingkey="xx_lda_prelim_score")
	fdr = CombinedAnalysisFDR(ndecoys=ndecoys,rankingkey="main_var_xx_swath_prelim_score")
        fdr.add_ids(csv.DictReader(open(target_file),dialect="excel-tab"))
    else:
        # fdr = SeparateAnalysisFDR(ndecoys=ndecoys,rankingkey="xx_lda_prelim_score")
        fdr = SeparateAnalysisFDR(ndecoys=ndecoys,rankingkey="main_var_xx_swath_prelim_score")
        fdr.add_target_ids(csv.DictReader(open(target_file),dialect="excel-tab"))
        for df in decoy_files:
	    fdr.add_decoy_ids(csv.DictReader(open(df),dialect="excel-tab"))

    for sc,qv,ta,de in fdr.allqvalues():
	print sc,qv,ta,de

    # for sc in range(10,-1,-1):
    #	print sc,fdr.qvalue(sc)

