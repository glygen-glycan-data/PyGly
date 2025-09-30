#!/bin/env python3.12

import os
import sys
import time
import glob
import findpygly
import pygly.alignment
from pygly.GlycanFormatter import GlycoCTFormat, WURCS20Format
from pygly.GlycanResource.GlyTouCan import GlyTouCanNoCache, GlyTouCan
from pygly.GlycanResource.GlyCosmos import GlyCosmosNoCache, GlyCosmos
from getwiki import GlycoMotifWiki
import alignmentindexchecker

import json
from collections import defaultdict
from pygly.GlycanFormatter import GlycoCTFormat, WURCS20Format, GlycanParseError
from operator import itemgetter

import argparse

from distproc import DistributedProcessing as dp

def check_idmaps(motifids,glycanids,idmaps): ###ID MAPS
    bad = 0
    for idmap in idmaps:
        idmapids = [ (t[0].id(),t[1].id()) for t in idmap ]
        bad = 0
        if set(map(itemgetter(0),idmapids)) != motifids:
            bad += 1
        if not set(map(itemgetter(1),idmapids)) <= glycanids:
            print(set(map(itemgetter(1),idmapids)),file=sys.stderr)
            print(glycanids,file=sys.stderr)
            print(set(map(itemgetter(1),idmapids)) <= glycanids,file=sys.stderr)
            bad += 2
        if len(set(map(itemgetter(1),idmapids))) != len(motifids):
            bad += 4
        if len(map(itemgetter(0),idmapids)) != len(set(map(itemgetter(0),idmapids))):
            bad += 8
        if len(map(itemgetter(1),idmapids)) != len(set(map(itemgetter(1),idmapids))):
            bad += 16
        if len(set(map(itemgetter(0),idmapids))) != len(set(idmapids)):
            bad += 32
        if bad > 0:
            print("Bad idmap:",bad,idmapids,file=sys.stderr)
            break
    return bad == 0

def idmaps_toids(idmaps,aligner):
    newidmaps = []
    for idmap in idmaps:
        idmapids = [ ti for t in idmap for ti in aligner.monoidmap(*t) ]
        newidmaps.append(idmapids)
    return newidmaps

def get_match_index (motifacc,structacc,idmaps,glycan=None):
    allstructids = set()
    allstructlinkids = set()
    for idmap in idmaps: 
        structids = set(t[1] for t in idmap)
        allstructids.update(structids)
        monoids = set(t1 for t1 in structids if '.' not in t1)
        linkids = set()
        if glycan:
            for l in glycan.all_links(uninstantiated=True):
                if l.parent().external_descriptor_id() in monoids and \
                   l.child().external_descriptor_id() in monoids:
                    if l.instantiated():
                        linkids.add((l.parent().external_descriptor_id(),l.child().external_descriptor_id()))
                    else:
                        linkids.add(('',l.child().external_descriptor_id()))
            if len(monoids) != (len(linkids) + 1):
                print("Warning: Bad linkids length motifacc: %s structacc: %s"%(motifacc,structacc),file=sys.stderr)
            allstructlinkids.update(linkids)
    x = "Y:"+",".join(str(i) for i in sorted(allstructids))
    if glycan:
        x += ":"+",".join("%s-%s"%p for p in sorted(allstructlinkids))
    return(x)

class Worker(object):

    def init(self,shared_data):
        self.wp = WURCS20Format()
        self.motif_gobjs = dict()
        for acc,seq in shared_data['motifs'].items():
            try:
                self.motif_gobjs[acc] = self.wp.toGlycan(seq)
            except GlycanParseError:
                continue
            except TypeError:
                continue
        self.task_blocks = list(shared_data['task_blocks'])

        self.nodes_cache = pygly.alignment.ConnectedNodesCache()

        self.loose_matcher = pygly.alignment.MotifInclusive(connected_nodes_cache=self.nodes_cache)
        self.loose_nred_matcher = pygly.alignment.NonReducingEndMotifInclusive(connected_nodes_cache=self.nodes_cache)
        
        self.strict_matcher = pygly.alignment.MotifStrict(connected_nodes_cache=self.nodes_cache)
        self.strict_nred_matcher = pygly.alignment.NonReducingEndMotifStrict(connected_nodes_cache=self.nodes_cache)

    def do_block_task(self,task,**kwargs):

        if not hasattr(self,'wp'):
            self.init(dict(kwargs.get('shared_data',{}).items()))
     
        lines = []
        for taski in self.task_blocks[task]['task_block']:
            lines.extend(self.do_task(taski,**kwargs))
        return lines
        
    def do_task(self,task,**kwargs):

        if not hasattr(self,'wp'):
            self.init(dict(kwargs.get('shared_data',{}).items()))
     
        glycan_acc = task['acc']
        glycan_seq = task['seq']
        try:
            glycan = self.wp.toGlycan(glycan_seq)
        except:
            return []

        self.nodes_cache.clear()

        start = time.time()
    
        lines = []
        for macc,motif in sorted(self.motif_gobjs.items()):

            idmaps_loose_core = []
            loose_core = self.loose_matcher.leq(motif, glycan, rootOnly=True, anywhereExceptRoot=False, underterminedLinkage=True, idmaps=idmaps_loose_core)
            idmaps_loose_core = idmaps_toids(idmaps_loose_core,self.loose_matcher)
        
            idmaps_loose_noncore = []
            loose_noncore = self.loose_matcher.leq(motif, glycan, rootOnly=False, anywhereExceptRoot=True, underterminedLinkage=True, idmaps=idmaps_loose_noncore)
            idmaps_loose_noncore = idmaps_toids(idmaps_loose_noncore,self.loose_matcher)
 
            loose_substructure = loose_core or loose_noncore
            idmaps_loose_substructure = list(idmaps_loose_core) + list(idmaps_loose_noncore)

            loose_whole = False
            idmaps_loose_whole = [] 
            if loose_core and self.loose_matcher.whole_glycan_match_check(motif, glycan):
                loose_whole=True
                idmaps_loose_whole = list(idmaps_loose_core)
            
            loose_nred=False
            idmaps_loose_nred = [] 
            if not motif.repeated() and not glycan.repeated() and loose_substructure:
                loose_nred = self.loose_nred_matcher.leq(motif, glycan, underterminedLinkage=True, idmaps=idmaps_loose_nred)
                idmaps_loose_nred = idmaps_toids(idmaps_loose_nred,self.loose_nred_matcher)

            strict_core, strict_noncore, strict_whole, strict_nred = False, False, False, False
            idmaps_strict_core = []
            idmaps_strict_noncore =[]
            idmaps_strict_whole = []
            idmaps_strict_nred = []
        
            if loose_core:
                strict_core = self.strict_matcher.leq(motif, glycan, rootOnly=True, anywhereExceptRoot=False, underterminedLinkage=False, idmaps=idmaps_strict_core)
                idmaps_strict_core = idmaps_toids(idmaps_strict_core,self.strict_matcher)
      
            if loose_noncore:
                strict_noncore = self.strict_matcher.leq(motif, glycan, rootOnly=False, anywhereExceptRoot=True, underterminedLinkage=False, idmaps=idmaps_strict_noncore)
                idmaps_strict_noncore = idmaps_toids(idmaps_strict_noncore,self.strict_matcher)
     
            strict_substructure = strict_core or strict_noncore
            idmaps_strict_substructure = list(idmaps_strict_core) + list(idmaps_strict_noncore) #partial changed here

            if strict_core and self.strict_matcher.whole_glycan_match_check(motif, glycan):
                strict_whole = True
                idmaps_strict_whole = list(idmaps_strict_core)
       
            if loose_nred and strict_substructure:
                strict_nred = self.strict_nred_matcher.leq(motif, glycan, underterminedLinkage=False, idmaps=idmaps_strict_nred)
                idmaps_strict_nred = idmaps_toids(idmaps_strict_nred,self.strict_nred_matcher)
            
            res0 = [loose_core, loose_substructure, loose_whole, loose_nred, 
                    strict_core, strict_substructure,strict_whole, strict_nred]
        
            if True not in res0: 
                continue

            index_numbers = [idmaps_loose_core,idmaps_loose_substructure,idmaps_loose_whole,idmaps_loose_nred, idmaps_strict_core, idmaps_strict_substructure , idmaps_strict_whole, idmaps_strict_nred]
        
            res1 = []
        
            for mt,idmaps in zip(res0,index_numbers):
                if mt:
                    indices=get_match_index(macc,glycan_acc,idmaps,glycan)
                    res1.append(indices)
                else:
                    res1.append("N")
        
            line = [macc, glycan_acc] + res1
            lines.append(line)

        # elapsed = time.time() - start
        # if elapsed >= 10:
        #     self.log("Analyze",glycan_acc,"elapsed: %.2f sec."%(elapsed,))

        return lines

if __name__ == "__main__":


    parser = argparse.ArgumentParser(description="Compute Alignments")
    parser.add_argument( '--smwenv', type = str, default = 'DEV', help = 'SMW instance. Default: DEV.')
    parser.add_argument( '-o', '--output', type = str, help = 'Output filename. Optional')
    parser.add_argument( '-A', '--accessions', type = str, nargs = "+", help = 'Motif accessions to align. Optional')
    parser.add_argument( '-q', '--quiet', action = 'store_true', default = False, help = 'No logging.')
    parser.add_argument( '-v', '--verbose', action = 'store_true', default = False, help = 'Verbose logging.')
    dp.add_arguments(parser)


    args = parser.parse_args()
    workers = dp.parse_args(parser)
    verbose = None
    if args.quiet:
        verbose = False
    if args.verbose:
        verbose = True

    worker = Worker()
    dp.start_if_worker(workers,worker.do_block_task)

    gtc = GlyTouCan(usecache=True)

    motifs = dict()
    if args.accessions is not None and len(args.accessions) > 0:
        for acc in args.accessions:
            motifs[acc] = gtc.getseq(acc,'wurcs')
    else:
        w = GlycoMotifWiki(smwenv=args.smwenv)
        for mid,acc,seq in w.itermotifgtcwurcs():
            if acc in motifs:
                continue
            motifs[acc]=seq

    archived = set()
    gco = GlyCosmosNoCache()
    for acc in gco.archived():
        archived.add(acc["accession"])
        
    gtcallseq = [ dict(seq=t[1],acc=t[0]) for t in sorted(filter(lambda t: t[0] not in archived,map(itemgetter(0,2),gtc.allseq(format="wurcs")))) ]

    blocksize = 100
    task_blocks = []
    for blockstart in range(0,len(gtcallseq),blocksize):
        task_blocks.append(dict(task_block=gtcallseq[blockstart:(blockstart+blocksize)]))
    tasks = list(range(0,len(task_blocks)))

    outhandle = None

    for res in dp.process(tasks=tasks,target=worker.do_block_task,workers=workers,verbose=verbose,motifs=motifs,task_blocks=task_blocks):
        if outhandle is None:
            outhandle = sys.stdout
            if args.output:
                outhandle = open(args.output,'w')
        for line in res:
            print("\t".join(line),file=outhandle)

    if args.output:
        outhandle.close()
