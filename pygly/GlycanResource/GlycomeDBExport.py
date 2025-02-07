
import sys, os, gzip

from .WebServiceResource import WebServiceResource

class GlycomeDBExport(WebServiceResource):
    apiurl = 'https://raw.githubusercontent.com/ReneRanzinger/org.glycomedb.export.glygen/refs/heads/main/export'
    def __init__(self,**kw):
        kw['iniFile'] = os.path.join(os.path.dirname(os.path.realpath(__file__)),"glycomedbexport.ini")
        super(GlycomeDBExport,self).__init__(**kw)

    def allgtc(self,resource=None):
        # "id","glytoucan","taxon","resource","resource_id"
        for row in self.query_taxmap():
            if resource == None or row['resource'] == resource:
                yield row['resource_id'],row['glytoucan'],row['resource']
            if resource == None or row['resource'] == 'glycomedb':
                yield row['id'],row['glytoucan'],'glycomedb'

    def alltaxa(self,resource=None):
         for row in self.query_taxmap():
            if resource == None or row['resource'] == resource:
               yield row['resource_id'],row['glytoucan'],row['taxon'],row['resource']

    def carbbank_allgtc(self):
        for row in self.query_taxmap():
            if row['resource'] == 'carbbank':
                try:
                    rsid = int(row['resource_id'])
                    yield row['resource_id'],row['glytoucan']
                except ValueError:
                    pass

    def carbbank_alltaxa(self):
        for row in self.query_taxmap():
            if row['resource'] == 'carbbank':
                try:
                    rsid = int(row['resource_id'])
                    yield rsid,row['glytoucan'],row['taxon']
                except ValueError:
                    pass

