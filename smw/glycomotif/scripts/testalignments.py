#!/bin/env python2
import os
import sys
import time
import findpygly
import pygly.alignment
from pygly.GlycanResource.GlyTouCan import GlyTouCanNoPrefetch

gtc = GlyTouCanNoPrefetch()

nodes_cache = pygly.alignment.ConnectedNodesCache()

loose_matcher = pygly.alignment.MotifInclusive(connected_nodes_cache=nodes_cache)
loose_nred_matcher = pygly.alignment.NonReducingEndMotifInclusive(connected_nodes_cache=nodes_cache)

strict_matcher = pygly.alignment.MotifStrict(connected_nodes_cache=nodes_cache)
strict_nred_matcher = pygly.alignment.NonReducingEndMotifStrict(connected_nodes_cache=nodes_cache)

motif_gobj = gtc.getGlycan(sys.argv[1],format='wurcs')
glycan_obj = gtc.getGlycan(sys.argv[2],format='wurcs')

# Loose match first
loose_core = loose_matcher.leq(motif_gobj, glycan_obj, rootOnly=True, anywhereExceptRoot=False, underterminedLinkage=True)
loose_substructure_partial = False
if not loose_core:
    loose_substructure_partial = loose_matcher.leq(motif_gobj, glycan_obj, rootOnly=False, anywhereExceptRoot=True, underterminedLinkage=True)
loose_substructure = loose_core or loose_substructure_partial

loose_whole = False
if loose_core and loose_matcher.whole_glycan_match_check(motif_gobj, glycan_obj):
    loose_whole = True

loose_nred = False
if not motif_gobj.repeated() and not glycan_obj.repeated() and loose_substructure:
    loose_nred = loose_nred_matcher.leq(motif_gobj, glycan_obj, underterminedLinkage=True)


# if inclusive, then try to match strict
strict_core, strict_substructure_partial, strict_whole, strict_nred = False, False, False, False
if loose_core:
    strict_core = strict_matcher.leq(motif_gobj, glycan_obj, rootOnly=True, anywhereExceptRoot=False, underterminedLinkage=False)

if loose_substructure_partial:
    if not strict_core:
        strict_substructure_partial = strict_matcher.leq(motif_gobj, glycan_obj, rootOnly=False, anywhereExceptRoot=True, underterminedLinkage=False)

if strict_core and strict_matcher.whole_glycan_match_check(motif_gobj, glycan_obj):
    strict_whole = True

strict_substructure = strict_core or strict_substructure_partial

if loose_nred and strict_substructure:
    strict_nred = strict_nred_matcher.leq(motif_gobj, glycan_obj, underterminedLinkage=False)

print "Motif:       ", sys.argv[1]
print "Glycan:      ", sys.argv[2]
print "Loose core:  ", loose_core
print "Loose subst: ", loose_substructure
print "Loose whole: ", loose_whole
print "Loose nred:  ", loose_nred
print "Strict core: ", strict_core
print "Strict subst:", strict_substructure
print "Strict whole:", strict_whole
print "Strict nred: ", strict_nred
