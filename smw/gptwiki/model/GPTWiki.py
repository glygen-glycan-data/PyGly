
__all__ = [ "GPTWiki", "Glycan", "Protein", "Peptide", "TransitionGroup", "Transition" ]

import sys
from operator import itemgetter

from smw import SMW

class Glycan(SMW.SMWClass):
    template = 'Glycan'
    classes = [ 'N-Glycan', 'O-Glycan' ]

    @staticmethod
    def pagename(**kwargs):
        assert kwargs.get('accession')
        return kwargs.get('accession')
    
    def asclass(self,value):
        if value not in self.classes:
	    raise ValueError("Bad class value: %r"%(value,))
	return value

    def toPython(self,data):
	data = super(Glycan,self).toPython(data)

        # name is newline separated
        if isinstance(data.get('name'),basestring):
            data['name'] = map(lambda s: s.strip(),data.get('name').split('\n'))

        # mw is a float
        if isinstance(data.get('mw'),basestring):
            data['mw'] = float(data.get('mw'))

        # class is comma separated, with specific possible values, sorted, so behaves as set
        if isinstance(data.get('class'),basestring):
            data['class'] = sorted(map(self.asclass,map(lambda s: s.strip(),data.get('class').split(','))))
	elif data.get('class') != None:
	    data['class'] = sorted(map(self.asclass,data.get('class')))
	
	return data

    def toTemplate(self,data):
	data = super(Glycan,self).toTemplate(data)

        if 'name' in data:
	    data['name'] = "\n".join(data['name'])

        if 'mw' in data:
            data['mw'] = str(data.get('mw'))
        
	if 'class' in data:
	    data['class'] = ",".join(sorted(data['class']))

	return data
	
class Protein(SMW.SMWClass):
    template = 'Protein'

    @staticmethod
    def pagename(**kwargs):
        assert kwargs.get('accession')
        return kwargs.get('accession')
    
class Transition(SMW.SMWClass):
    template = 'Transition'

    def toPython(self,data):
	data = super(Transition,self).toPython(data)

        for k in ('nrt','rt','mz1','mz2'):
            if isinstance(data.get(k),basestring):
                data[k] = float(data.get(k))
                
        for k in ('z1','z2'):
            if isinstance(data.get(k),basestring):
                data[k] = int(data.get(k))

	return data

    def toTemplate(self,data):
	data = super(Transition,self).toTemplate(data)

        for k in ('nrt','rt','mz1','mz2'):
            if k in data:
                data[k] = "%.3f"%(data[k])

        for k in ('z1','z2'):
            if k in data:
                data[k] = "%d"%(data[k])

	return data

class TransitionGroup(SMW.SMWClass):
    template = 'TransitionGroup'

    def astrans(self,tr):
        str = tr.split(';')
        return (str[0],float(str[1]))

    def asscans(self,sc):
        str = sc.split(';')
	if len(str) == 2:
            return (int(str[0]),str[1])
	if len(str) == 4:
            return (int(str[0]),str[1],float(str[2]),int(str[3]))
        return (int(str[0]),str[1],float(str[2]),int(str[3]),float(str[4]))

    def toPython(self,data):
	data = super(TransitionGroup,self).toPython(data)

        # name is newline separated
        if isinstance(data.get('transitions'),basestring):
            data['transitions'] = map(lambda t: self.astrans(t),data.get('transitions').split(','))

        if isinstance(data.get('scans'),basestring):
            data['scans'] = map(lambda t: self.asscans(t),data.get('scans').split(','))

        for k in ('rt','nrt','mz1'):
            if isinstance(data.get(k),basestring):
                data[k] = float(data.get(k))
                
        for k in ('z1','ntransition'):
            if isinstance(data.get(k),basestring):
                data[k] = int(data.get(k))

	return data

    @staticmethod
    def scan2str(scan):
	print scan
	if len(scan) == 2:
	    return "%d;%s"%(int(scan[0]),scan[1])
	if len(scan) == 4:
            return "%d;%s;%.3f;%d"%(int(scan[0]),scan[1],float(scan[2]),int(scan[3]))
        return "%d;%s;%.3f;%d;%.3f"%(int(scan[0]),scan[1],float(scan[2]),int(scan[3]),float(scan[4]))

    def toTemplate(self,data):
	data = super(TransitionGroup,self).toTemplate(data)

        if 'transitions' in data:
            data['transitions'] = ",".join(map(lambda t: "%s;%.1f"%t,data['transitions']))

        if 'scans' in data:
            data['scans'] = ",".join(map(self.scan2str,data['scans']))

        for k in ('rt','nrt','mz1'):
            if k in data:
                data[k] = "%.3f"%(data[k])

        for k in ('z1','ntransition'):
            if k in data:
                data[k] = "%d"%(data[k])

	return data


class Peptide(SMW.SMWClass):
    template = 'Peptide'

    @staticmethod
    def key(sequence,glycan,mod):
        l = [ sequence ]
	for t in sorted(glycan,key=lambda t: int(t[1][1:])):
	    l.append(":".join(map(str,t)))
	for t in sorted(mod,key=lambda t: (int(t[1][1:]),"%+.3f"%t[0])):
	    l.append("%+.3f"%(t[0],)+":"+":".join(t[1:]))
  	return "^".join(l)

    def asglycan(self,g):
        sg = g.split(';')
        return tuple(sg)

    def asmod(self,m):
        sm = m.split(';')
        return tuple([float(sm[0])] + sm[1:])

    def asalign(self,a):
        sa = a.split(';')
        return (sa[0],int(sa[1]),int(sa[2]))

    def toPython(self,data):
	data = super(Peptide,self).toPython(data)

        # mw is a float
        if isinstance(data.get('mw'),basestring):
            data['mw'] = float(data.get('mw'))
            
        if isinstance(data.get('glycan'),basestring):
            data['glycan'] = map(self.asglycan,data.get('glycan').split(','))

        if isinstance(data.get('mod'),basestring):
            data['mod'] = map(self.asmod,data.get('mod').split(','))

        if isinstance(data.get('alignment'),basestring):
            data['alignment'] = sorted(map(self.asalign,data.get('alignment').split(',')))
            # data['protein'] = sorted(set(map(itemgetter(0),data.get('alignment'))))

        if isinstance(data.get('transgroup'),basestring):
            data['transgroup'] = sorted(data.get('transgroup').split(','))

        if isinstance(data.get('sample'),basestring):
            data['sample'] = sorted(map(string.strip,data.get('sample').split('\n')))

	return data

    def toTemplate(self,data):
	data = super(Peptide,self).toTemplate(data)

        if 'mw' in data:
            data['mw'] = "%.3f"%(data.get('mw'),)

        if 'glycan' in data:
            data['glycan'] = ",".join(map(lambda t: "%s;%s"%t,sorted(data['glycan'],key=lambda t: int(t[1][1:]))))            

        if 'mod' in data:
            data['mod'] = ",".join(map(lambda t: "%+.3f;%s"%(t[0],";".join(t[1:])),
                                       sorted(data['mod'],key=lambda t: int(t[1][1:]))))

        if 'alignment' in data:
            data['alignment'] = ",".join(map(lambda t: "%s;%d;%d"%t,sorted(data['alignment'])))
            # data['protein'] = ",".join(set(map(lambda t: t.split(';')[0],sorted(data['alignment'].split(',')))))

        if 'transgroup' in data:
            data['transgroup'] = ",".join(data.get('transgroup'))

        if 'sample' in data:
            data['sample'] = "\n".join(data.get('sample'))

	return data    
	
class GPTWiki(SMW.SMWSite):
    _name = 'gptwiki'

    def __init__(self):
        super(GPTWiki,self).__init__()
        self._peptides = {}
        for p in map(self.get,self.itercat('Peptide')):
            key = Peptide.key(p.get('sequence'),p.get('glycan',[]),p.get('mod',[]))
            self._peptides[key] = p.get('id')
        self._transgroups = {}
        for tg in map(self.get,self.itercat('TransitionGroup')):
            key = (tg.get('peptide'),tg.get('z1'),tg.get('spectra'),tg.get('method'))
            self._transgroups[key] = tg.get('id')
        self._transitions = {}
        for t in map(self.get,self.itercat('Transition')):
            key = (t.get('peptide'),t.get('z1'),t.get('label'),t.get('z2'))
            self._transitions[key] = t.get('id')

    def addprotein(self,accession,**kw):
        p = self.get(accession)
        if p:
            p.update(kw)
        else:
            p = Protein(accession=accession,**kw)
        self.put(p)
        
    def addglycan(self,accession,**kw):
        p = self.get(accession)
        if p:
            p.update(kw)
        else:
            p = Glycan(accession=accession,**kw)
        self.put(p)

    def newtransid(self):
        maxid = 0
        for i,id in enumerate(sorted(map(lambda id: int(id[2:]),self._transitions.values()))):
            if (i+1) < id:
                return "TR%06d"%(i+1,)
            maxid = id
        return "TR%06d"%(maxid+1,)

    def addtransition(self,peptide,z1,label,z2,**kw):
        id = self._transitions.get((peptide,z1,label,z2))
        if id:
            t = self.get(id)
            t.update(**kw)
        else:
            id = self.newtransid()
            t = Transition(id=id,peptide=peptide,z1=z1,label=label,z2=z2,**kw)
            self._transitions[(peptide,z1,label,z2)] = id
        return t,self.put(t)

    def newtransgrpid(self):
        maxid = 0
        for i,id in enumerate(sorted(map(lambda id: int(id[2:]),self._transgroups.values()))):
            if (i+1) < id:
                return "TG%06d"%(i+1,)
            maxid = id
        return "TG%06d"%(maxid+1,)

    def addtransgroup(self,peptide,z1,spectra,method,**kw):
        id = self._transgroups.get((peptide,z1,spectra,method))
        # print self._transgroups
        if id:
            tg = self.get(id)
            tg.update(**kw)
        else:
            id = self.newtransgrpid()
            tg = TransitionGroup(id=id,peptide=peptide,z1=z1,spectra=spectra,method=method,**kw)
            self._transgroups[(peptide,z1,spectra,method)] = id
        return tg,self.put(tg)

    def iterpeptides(self):
	for p in map(self.get,self.itercat('Peptide')):
            yield p

    def newpeptideid(self):
        maxid = 0
        for i,id in enumerate(sorted(map(lambda id: int(id[2:]),self._peptides.values()))):
            if (i+1) < id:
                return "PE%06d"%(i+1,)
            maxid = id
        return "PE%06d"%(maxid+1,)

    def findpeptide(self,sequence,glycans=[],mods=[]):
        key = Peptide.key(sequence,glycan=glycans,mod=mods)
        return key,self._peptides.get(key)

    def addpeptide(self,sequence,glycans=[],mods=[],**kw):
        key,id = self.findpeptide(sequence,glycans=glycans,mods=mods)
        if id:
            p = self.get(id)
            p.update(**kw)
        else:
            id = self.newpeptideid()
            p = Peptide(id=id,sequence=sequence,glycan=glycans,mod=mods,**kw)
            self._peptides[key] = id
	return p,self.put(p)
            

