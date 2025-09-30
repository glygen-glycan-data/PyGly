#!/bin/env python3
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
loose_core_idmaps = []
loose_core = loose_matcher.leq(motif_gobj, glycan_obj, rootOnly=True, anywhereExceptRoot=False, underterminedLinkage=True, idmaps=loose_core_idmaps)
print(loose_core, loose_core_idmaps)
loose_substructure_partial_idmaps = []
loose_substructure_partial = loose_matcher.leq(motif_gobj, glycan_obj, rootOnly=False, anywhereExceptRoot=True, underterminedLinkage=True,idmaps=loose_substructure_partial_idmaps)
loose_substructure = (loose_core or loose_substructure_partial)
loose_substructure_idmaps = (loose_core_idmaps + loose_substructure_partial_idmaps)

loose_whole = False
loose_whole_idmaps = []
if loose_core and loose_matcher.whole_glycan_match_check(motif_gobj, glycan_obj):
    loose_whole = True
    loose_whole_idmaps = loose_core_idmaps

loose_nred = False
loose_nred_idmaps = []
if not motif_gobj.repeated() and not glycan_obj.repeated() and loose_substructure:
    loose_nred = loose_nred_matcher.leq(motif_gobj, glycan_obj, underterminedLinkage=True, idmaps=loose_nred_idmaps)

# if inclusive, then try to match strict
strict_core, strict_substructure, strict_substructure_partial, strict_whole, strict_nred = False, False, False, False, False
strict_core_idmaps, strict_substructure_idmaps, strict_substructure_partial_idmaps, strict_whole_idmaps, strict_nred_idmaps = [], [], [], [], []
if loose_core:
    strict_core = strict_matcher.leq(motif_gobj, glycan_obj, rootOnly=True, anywhereExceptRoot=False, underterminedLinkage=False,idmaps=strict_core_idmaps)

if loose_substructure_partial:
    strict_substructure_partial = strict_matcher.leq(motif_gobj, glycan_obj, rootOnly=False, anywhereExceptRoot=True, underterminedLinkage=False,idmaps=strict_substructure_partial_idmaps)

strict_substructure = (strict_core or strict_substructure_partial)
strict_substructure_idmaps = (strict_core_idmaps + strict_substructure_partial_idmaps)

if strict_core and strict_matcher.whole_glycan_match_check(motif_gobj, glycan_obj):
    strict_whole = True
    strict_whole_idmaps = strict_core_idmaps

if loose_nred and strict_substructure:
    strict_nred = strict_nred_matcher.leq(motif_gobj, glycan_obj, underterminedLinkage=False,idmaps=strict_nred_idmaps)

def tostridset(idmaps):
    ids = set()
    for idmap in idmaps:
        ids.update([p[1].id() for p in idmap])
    return ",".join(map(str,sorted(ids)))

print("Motif:               ", sys.argv[1])
print("Glycan:              ", sys.argv[2])
print("Loose core:          ", loose_core, tostridset(loose_core_idmaps))
print("Loose subst partial: ", loose_substructure_partial, tostridset(loose_substructure_partial_idmaps))
print("Loose subst:         ", loose_substructure, tostridset(loose_substructure_idmaps))
print("Loose nred:          ", loose_nred, tostridset(loose_nred_idmaps))
print("Loose whole:         ", loose_whole, tostridset(loose_whole_idmaps))
print("Strict core:         ", strict_core, tostridset(strict_core_idmaps))
print("Strict subst partial:", strict_substructure_partial, tostridset(strict_substructure_partial_idmaps))
print("Strict subst:        ", strict_substructure, tostridset(strict_substructure_idmaps))
print("Strict whole:        ", strict_whole, tostridset(strict_whole_idmaps))
print("Strict nred:         ", strict_nred, tostridset(strict_nred_idmaps))
