
import sys
from ReferenceTable import ReferenceTable
from Monosaccharide import *


class MonoFactory(ReferenceTable):
    def new(self,key):
        return self[key].clone()
    def parseSection(self,name,kv):
	try:
            m = Monosaccharide()
            m.set_id(name)
            aliases = [name]
            try:
                cls,const = constantLookup(kv['anomer'])
                m.set_anomer(const)
            except KeyError:
                pass
            cls,const = constantLookup(kv['superclass'])
            m.set_superclass(const)
            try:
                m.set_ring_start(int(kv['ring_start']))
            except KeyError:
                pass
            try:
                m.set_ring_end(int(kv['ring_end']))
            except KeyError:
                pass
            m.set_stem(*map(lambda st: constantLookup(st)[1],kv['stem'].split()))
            try:
                m.set_config(*map(lambda cg: constantLookup(cg)[1],kv['config'].split()))
            except KeyError:
                pass
            mods = kv.get('mods','').split()
            for i in range(0,len(mods),2):
                m.add_mod(mods[i],constantLookup(mods[i+1])[1])
            substs = filter(None,map(str.strip,kv.get('substituent','').split(';')))
            for subst in substs:
                split_subst = subst.split()
                const = constantLookup(split_subst[0])[1]
                kwargs = {}
                for i in range(1,len(split_subst),2):
                    if split_subst[i].endswith('_pos'):
                        kwargs[split_subst[i]] = int(split_subst[i+1])
                    else:
                        kwargs[split_subst[i]] = constantLookup(split_subst[i+1])[1]
                m.add_substituent(const,**kwargs)
            aliases.extend(map(str.strip,kv.get('aliases','').split(';')))
            aliases = filter(None,aliases)
	except:
	    print >>sys.stderr, "Problem with section", name
	    raise
        return [(a,m) for a in aliases]

        
