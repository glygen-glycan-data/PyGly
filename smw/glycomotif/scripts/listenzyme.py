#!/bin/env python27

from getwiki import GlycoMotifWiki

w = GlycoMotifWiki()

for e in w.iterenzyme():
    # if e.has('dscription'):
    #     e.delete('dscription')
    # if e.has('disease_id'):
    #     e.delete('disease_id')
    # if e.has('phenotype_soure_key'):
    #     e.set('phenotype_source_key',e.get('phenotype_soure_key'))
    #     e.delete('phenotype_soure_key')
    if 'disease' in e.get('description'):
        desc = e.get('description').split('\n')[0]
        e.set('description',desc)
    if w.put(e):
        print e.get('genename')
