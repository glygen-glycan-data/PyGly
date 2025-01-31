
from __future__ import print_function

from . Monosaccharide import Anomer, Stem, SuperClass, Config, Substituent, Mod

import re, copy

class Manipulation(object):

    def __call__(self,g):
        return self.make(g)

    def make(self,g):
        g1 = g.clone()
        self.manip(g1)
        return g1

class RemoveAlditol(Manipulation):

    def manip(self,g):
        assert g.has_root()
        r = g.root()
        r.remove_mod(Mod.aldi,(1,))

class Archetype(Manipulation):
 
    def manip(self,g):
        # fully_defined does not require the root monosaccharide anomer
        # or ring start and end be defined
        assert(g.has_root() and not g.repeated())
        r = g.root()
        r.set_anomer(Anomer.missing)
        r.set_ring_start(None)
        r.set_ring_end(None)
        r.remove_mod(Mod.aldi,(1,))

class Topology(Manipulation):

    def manip(self,g):
        for m in g.all_nodes():
            for l in m.links_with_uninstantiated():

                lparent = l.parent()  # might be monosaccharide or substituent

                if lparent.is_monosaccharide():
                    if l.parent_pos() != set([m.ring_start()]):
                        l.set_parent_pos(None)
                else:
                    # Substituent in link
                    upper_link = l.parent().any_parent_link()
                    #if upper_link.parent_pos() != set([m.ring_start()]):
                    upper_link.set_parent_pos(None)

                    if l.child_pos() != set([m.ring_start()]):
                        l.set_child_pos(None)
                # if l.parent_pos() != set([m.ring_start()]):
                #    l.set_parent_pos(None)
                # if l.parent_pos() == set([1]) and m.ring_start() != 1 and l.child_pos() == l.child().ring_start():
                #     l.set_parent_pos(None)
                # if l.parent_pos() != set([1]) and l.parent_pos() != set([m.ring_start()]):
                #     l.set_parent_pos(None)
                # if l.child_pos() != set([1]) and l.child().superclass() in (SuperClass.HEX,):
                #     l.set_child_pos(None)
                if not l.instantiated() and l.child_pos() == set([1]) and l.child().superclass() in (SuperClass.NON,SuperClass.OCT):
                    l.set_child_pos(None)
            if m.anomer() != Anomer.uncyclized:
                m.set_anomer(Anomer.missing)
        undets = list(g.undetermined_roots())
        g.set_undetermined(undets)

class Composition(Manipulation):

    topo = Topology()

    non_floating_substs = set(map(lambda k: getattr(Substituent,k),"""
        amino
        nAcetyl
        anhydro
        fluoro
        nformyl
        chloro
        ndimethyl
        hydroxymethyl
        nsulfate
        bromo
        nglycolyl
        methyl_oxygen_lost
        scarboxyethyl
        thio
        namidino
        iodo
        nsuccinate
        acetyl_oxygen_lost
        ethanolamine
        nmethyl
        sulfate_oxygen_lost
        spyruvate
        phosphate_oxygen_lost
    """.split()))
    floating_substs = set(map(lambda k: getattr(Substituent,k),"""
        methyl
        acetyl
        formyl
        sulfate
        phosphate
        slactate
        rlactate
        rcarboxyethyl
        phosphoethanolamine
        pyrophosphate
        phosphocholine
        diphosphoethanolamine
        ethyl
        amino_oxygen_preserved
        glycolyl
        xlactate
        triphosphate
        phosphate_bridged
        pyruvate
    """.split()))

    def manip(self,g):
        self.topo.manip(g)
        undets = []
        for m in g.all_nodes(undet_subst=True):
            undets.append(m)
        floating = []
        for m in undets:
            m.set_connected(False)
            m.clear_links()
            m.clear_parent_links()
            if m.is_monosaccharide():
                if "%s:%s"%(m.ring_start(),m.ring_end()) == "0:0" and m.anomer() == Anomer.uncyclized:
                    pass
                else:
                    m.set_ring_start(None)
                    m.set_ring_end(None)
            for sublink in list(m.substituent_links()):
                assert sublink.child().name() in self.floating_substs or sublink.child().name() in self.non_floating_substs, str(sublink.child())
                if sublink.child().name() in self.floating_substs:
                    m.remove_substituent_link(sublink)
                    sublink.child().del_parent_link(sublink)
                    sublink.child().set_connected(False)
                    floating.append(sublink.child())
                    sublink.child().clear_links()
        undets.extend(floating)
        if len(undets) == 1:
            g.set_root(undets[0])
            g.set_undetermined([])
            g.root().set_connected(True)
        else:
            g.set_root(None)
            g.set_undetermined(undets)

        # for m in g.all_nodes():
        #     print m, m._links, map(str,m.parent_links())

class BaseComposition(Manipulation):

    comp = Composition()

    def manip(self,g):
        self.comp.manip(g)
        for m in g.all_nodes():
            m.set_stem(Stem.missing)
            m.set_config(Config.missing)
            m.set_ring_start(None)
            m.set_ring_end(None)
            m.set_anomer(Anomer.missing)

class LevelSniffer(object):
    topology = Topology()
    composition = Composition()
    basecomposition = BaseComposition()
    
    def __call__(self,g):
        return self.sniff(g)

    def sniff(self,g):
        if g.has_root():
            if g.equals(self.topology(g)):
                return "Topology"
            else:
                return "Saccharide"
        else:
            if g.equals(self.basecomposition(g)):
                return "BaseComposition"
            elif g.equals(self.composition(g)):
                return "Composition"
        raise ValueError

class WURCSManipulation(object):
    def __init__(self):
        self.alpha2ind = dict()
        self.ind2alpha = dict()
        for i in range(ord('a'),ord('z')+1):
            self.alpha2ind[chr(i)] = i-ord('a')+1
            self.ind2alpha[i-ord('a')+1] = chr(i)
        for i in range(ord('A'),ord('Z')+1):
            self.alpha2ind[chr(i)] = i-ord('A')+1+26
            self.ind2alpha[i-ord('A')+1+26] = chr(i)

    def __call__(self,seq):
        return self.modify(seq)

    @staticmethod
    def intorstr(s):
        try:
            return int(s)
        except ValueError:
            pass
        return s

    def deconstruct(self,seq):
        head,cnts,rest = seq.split('/',2)
        monos,rest = rest.split(']/',1)
        inds,links = rest.split('/',1)
                                                                                                           
        cnts = list(map(self.intorstr,cnts.split(',')))
        inds = list(map(int,inds.split('-')))
        monos = monos.lstrip('[').rstrip(']').split('][')
        monos = [ self.monodeconstruct(m) for m in monos ]
        links = filter(None,links.split('_'))
        links = [ self.linkdeconstruct(l) for l in links ]
        return dict(head=head,counts=cnts,monos=monos,indices=inds,links=links)

    def toseq(self,parts):
        head = parts['head']
        monos = parts['monos']
        inds = parts['indices']
        links = parts['links']
        monos,inds = self.repairindices(monos,inds)

        cnts = ",".join(map(str,[len(monos),len(inds),len(links)]))
        monos = '[' + ']['.join([ self.monoseq(m) for m in monos]) + ']'
        inds = '-'.join(map(str,inds))
        links = '_'.join([ self.linkseq(l) for l in links])
        return '/'.join([head,cnts,monos,inds,links])

    def monodeconstruct(self,monoseq):
        sm = monoseq.split('_')
        m1 = re.search(r'-\d[abxud]$',sm[0])
        if m1:
            base,anomer = sm[0].split('-')
            childpos,anomer = tuple(anomer)
            childpos = int(childpos)
            if len(sm) > 1:
                m2 = re.search(r'^\d-[0-9?]$',sm[1])
                if m2:
                    ring = sm[1]
                    substs=sm[2:]
                else:
                    ring = None
                    substs=sm[1:]
            else:
                ring = None
                substs=sm[1:]
        else:
            base = sm[0]
            anomer = None
            childpos = None
            ring = None
            substs=sm[1:]
        substs1 = []
        for s in substs:
            try:
                pos,sub = s.split('*',1)
            except ValueError:
                pos = s; sub = ""
            substs1.append(dict(pos=self.intorstr(pos),subst=sub))
        substs = substs1 
        return dict(skeleton=base,anomer=anomer,childpos=childpos,ring=ring,substituents=substs,wurcs=monoseq)

    def monoseq(self,m):
        mseq = m['skeleton']
        if m.get('anomer') and m.get('childpos'):
            mseq += '-' + str(m['childpos']) + m['anomer']
        if m.get('ring'):
            mseq += "_" + m['ring']
        if len(m['substituents']) > 0:
            mseq += "_" + "_".join(map(self.substseq,(m['substituents'])))
        if 'wurcs' in m:
            assert mseq == m['wurcs']
        return mseq

    def substseq(self,s):
        if s.get('subst'):
            return "%s*%s"%(s.get('pos'),s.get('subst'))
        return "%s"%(s.get('pos'),)

    linkre0 = re.compile(r'^(([a-zA-Z])(\d|\?))-(([a-zA-Z])(\d|\?))(~(n|(\d+)-(\d+)))?$')
    linkre1 = re.compile(r'^(([a-zA-Z])(\d|\?)(\|(([a-zA-Z])(\d|\?)))+)\}-\{(([a-zA-Z])(\d|\?)(\|(([a-zA-Z])(\d|\?)))+)$')
    linkre2 = re.compile(r'^(([a-zA-Z])(\d|\?)(\|(([a-zA-Z])(\d|\?)))+)\}\*(.*)$')
    linkre3 = re.compile(r'^(([a-zA-Z])(\d|\?))-(([a-zA-Z])(\d|\?))\*(.*?)(~(n|(\d+)-(\d+)))?$')
    def linkdeconstruct(self,linkseq):
        m = self.linkre0.search(linkseq)
        if m:
            parents = [ self.linkmonoparse(m.group(1)) ]
            children = [ self.linkmonoparse(m.group(4)) ]
            return dict(parents=parents,children=children,repeatlink=m.group(8))

        m = self.linkre1.search(linkseq)
        if m:
            parents = list(map(self.linkmonoparse,m.group(1).split('|')))
            children = list(map(self.linkmonoparse,m.group(8).split('|')))
            return dict(parents=parents,children=children)
            
        m = self.linkre2.search(linkseq)
        if m:
            parents = list(map(self.linkmonoparse,m.group(1).split('|')))
            subst=m.group(8)
            return dict(parents=parents,subst=subst)
            
        m = self.linkre3.search(linkseq)
        if m:
            parents = [ self.linkmonoparse(m.group(1)) ]
            children = [ self.linkmonoparse(m.group(4)) ]
            subst=m.group(7)
            repeatlink=m.group(9)
            return dict(parents=parents,children=children,subst=subst,repeatlink=repeatlink)
            
        return dict(seq=linkseq)

    def linkmonoparse(self,linkmono):
        try:
            pos = int(linkmono[1])    
        except ValueError:
            pos = linkmono[1]
        index = self.alpha2ind[linkmono[0]]
        return dict(pos=pos,index=index)

    def linkseq(self,l):
        if 'seq' in l:
            return l['seq']
        return "%s%s-%s%s"%(self.ind2alpha[l['parentindex']],l['parentpos'],
                            self.ind2alpha[l['childindex']],l['childpos'])

    def repairindices(self,monos,indices):
        unused = [ i for i in range(len(monos)) if (i+1) not in indices ]
        monoseq = [ self.monoseq(m) for m in monos ]
        newmonos = list(monos)
        newmonoseq = [ self.monoseq(m) for m in newmonos ]
        for i in range(len(newmonos)-1,-1,-1):
            try:
                ind = newmonoseq[:i].index(newmonoseq[i])
            except ValueError:
                if i not in unused:
                    continue
            del newmonos[i]
            del newmonoseq[i]
        indsmap = dict()
        for i in range(len(monos)):
            if i not in unused:
                indsmap[i+1] = newmonoseq.index(monoseq[i])+1
        inds = [ indsmap[ind] for ind in indices ]
        indsremap = dict()
        indsinvremap = dict()
        for j,i in enumerate(sorted(list(range(1,len(newmonoseq)+1)),key=lambda i: inds.index(i))):
            indsremap[j+1] = i
            indsinvremap[i] = (j+1)
        newmonoseq1 = []
        newmonos1 = []
        for i in range(len(newmonoseq)):
            newmonoseq1.append(newmonoseq[indsremap[i+1]-1])
            newmonos1.append(newmonos[indsremap[i+1]-1])
        inds1 = []
        for i in range(len(inds)):
            inds1.append(indsinvremap[inds[i]])
        newmonoseq = newmonoseq1
        newmonos = newmonos1
        inds = inds1
        assert(len(set(inds)) == len(newmonoseq) and min(inds) == 1 and max(inds) == len(newmonoseq))
        return newmonos,inds

class WURCSArchetype(WURCSManipulation):

    def mapskel(self,skel):
        m = re.search(r'^([hA][1234dx])[aUO]([1234dx]+[mhA])$',skel)
        if m:
            return m.group(1)+'U'+m.group(2)
        m = re.search(r'^([hA])[aUO]([1234dx]+[mhA])$',skel)
        if m:
            return m.group(1)+'U'+m.group(2)
        m = re.search(r'^[ahou]([1234568defxzEFZ]+[mhA])$',skel)
        if m:
            return 'u'+m.group(1)
        return None

    def redend_mono(self,monostr):
        redendmono = self.monodeconstruct(monostr)
        redendmono['skeleton'] = self.mapskel(redendmono['skeleton'])
        redendmono['anomer'] = None
        redendmono['childpos'] = None
        redendmono['ring'] = None
        if 'wurcs' in redendmono:
            del redendmono['wurcs']
        return self.monoseq(redendmono)

    def modify(self,seq):
        parts = self.deconstruct(seq)
        redendmono = parts['monos'][0]
        redendmono = copy.deepcopy(redendmono)
        skel2 = self.mapskel(redendmono['skeleton'])
        assert skel2, redendmono['skeleton']
        # assert redendmono['skeleton'] in self.skelmap, redendmono['skeleton']
        # skel1 = self.skelmap.get(redendmono['skeleton'])
        # assert skel1 == skel2, " ".join(map(str,[redendmono['skeleton'],skel1, skel2]))
        redendmono['skeleton'] = skel2
        redendmono['anomer'] = None
        redendmono['childpos'] = None
        redendmono['ring'] = None
        parts['monos'].insert(0,redendmono)
        parts['indices'] = [ i+1 for i in parts['indices'] ]
        parts['indices'][0] = 1
        return self.toseq(parts)

def showdiff(g,g2):
    gctstr = glycoct_parser.toStr(g)
    gctstr2 = glycoct_parser.toStr(g2)

    for l1,l2 in zip(gctstr.split('\n'),gctstr2.split('\n')):
        char = " "
        if l1 != l2:
            char = "!"
        print("%-40s"%(l1[:40],),char,"%-40s"%(l2[:40],),file=sys.stderr)

    print(sum(1 for _ in g.all_nodes(subst=True)))
    print(sum(1 for _ in g2.all_nodes(subst=True)))
    for m1,m2 in zip(g.all_nodes(subst=True),g2.all_nodes(subst=True)):
        print(str(m1))
        print(str(m2))

if __name__ == "__main__":

    import sys, csv, zipfile, traceback

    from GlycanFormatter import WURCS20Format, GlycoCTFormat, WURCS20ParseError, UnsupportedLinkError,  \
                                CircularError, GlycoCTParseError, MonoOrderLinkError, GlycanParseError
    from WURCS20MonoFormatter import UnsupportedMonoError

    # to run this testing framework, use the following command(s) for
    # human, mouse, and all GlyTouCan Saccharides.

    # setenv GLYTOUCAN "~/projects/GlyGen/GlyTouCan/current"
    #
    # egrep -w '(GlyTouCanAccession|Saccharide)' $GLYTOUCAN/comp.txt | \
    #     fgrep -f $GLYTOUCAN/humanbycomp.txt | \
    #     python27 manipulation.py $GLYTOUCAN/wurcs.zip $GLYTOUCAN/glycoct.zip 
    #
    # egrep -w '(GlyTouCanAccession|Saccharide)' $GLYTOUCAN/comp.txt | \
    #     fgrep -f $GLYTOUCAN/mousebycomp.txt | \
    #     python27 manipulation.py $GLYTOUCAN/wurcs.zip $GLYTOUCAN/glycoct.zip 
    #
    # egrep -w '(GlyTouCanAccession|Saccharide)' $GLYTOUCAN/comp.txt | \
    #     python27 manipulation.py $GLYTOUCAN/wurcs.zip $GLYTOUCAN/glycoct.zip 
    #

    wurcs_parser = WURCS20Format()
    glycoct_parser = GlycoCTFormat()

    topology = Topology()
    composition = Composition()
    basecomposition = BaseComposition()
    level = LevelSniffer()
    
    zf = zipfile.ZipFile(sys.argv[1])
    zf1 = zipfile.ZipFile(sys.argv[2])
    
    reader = csv.DictReader(sys.stdin,dialect='excel-tab')

    j = 0
    for i,d in enumerate(reader):
        acc = d['GlyTouCanAccession']
        typ = d['GlyTouCanType']
        topo = d.get('Topology')
        comp = d.get('Composition')
        bcomp = d.get('BaseComposition')

        print("\t".join(map(str,[i,j,acc,typ,topo,comp,bcomp])),file=sys.stderr)

        try:
            gtopo = None; gcomp = None; 
            g = wurcs_parser.toGlycan(zf.read(acc+'.txt'))
            if topo:
                gtopo = wurcs_parser.toGlycan(zf.read(topo+'.txt'))
            if comp:
                gcomp = wurcs_parser.toGlycan(zf.read(comp+'.txt'))
        except KeyError:
            print("Skip: Accession not in zip file",file=sys.stderr)
            # accession not in zip file
            continue
        except UnsupportedLinkError:
            print("Skip: Unsupported link error",file=sys.stderr)
            # unexpected link structure
            continue
        except MonoOrderLinkError:
            print("Skip: Unsupported monosaccharide order",file=sys.stderr)
            # link order issue
            continue
        except CircularError:
            print("Skip: Unsupported circular glycan",file=sys.stderr)
            # dodge this case...
            continue
        except UnsupportedMonoError as e:
            print("Skip:",e.message,file=sys.stderr)
            continue
        except GlycanParseError as e:
            print("Skip:",e.message,file=sys.stderr)
            continue


        if not gtopo or not gcomp:
            continue

        try:
            g3 = glycoct_parser.toGlycan(zf1.read(acc+'.txt'))
        except (KeyError, GlycoCTParseError):
            g3 = None

        if acc in ('G24172ZD',):
            # this is a (1+1) two monosaccharide case in which the
            # GlyTouCan topology has (-1+1).  Other (1+1) two
            # monosaccharide caes have GlyTouCan topology with (1+1) -
            # and this second decision seems to be more common.
            # Can't figure out how to decide whicoh one to do.
            print("Skip: (1+1) two monosaccharide issue",file=sys.stderr)
            continue

        if False and acc in ('G64632PP',):
            # the toplogy G83908ZR breaks my equality algorithm when compared to the 
            # topology generated from G64632PP - multiple valid id
            # mappings are possible, but the undetermined subtree only
            # captures one of them. Do I need to enumerate id mappings? 
            continue

        if False and acc in ('G00526WX','G03352IE','G04627JF','G04923UJ','G06441HU','G06602IJ','G07875MT','G08549SV','G08864PG',
                   'G10665RH','G10984JH','G12178PQ','G12612TV'):
            # WURCS parser problem w/ 2-5, 3-6 substituent. 
            continue

        print(glycoct_parser.toStr(g),file=sys.stderr)

        g1 = g.clone()

        print(glycoct_parser.toStr(g1),file=sys.stderr)

        assert g.equals(g1)
        assert g1.equals(g)

        gctstr = glycoct_parser.toStr(g)

        g2 = glycoct_parser.toGlycan(gctstr)

        print(glycoct_parser.toStr(g2),file=sys.stderr)

        # assert g.equals(g2)

        # assert g1.equals(g2)

        if acc not in ('G00727XI','G08545ST','G09485LV','G11795CV','G13902SG','G14164XT','G17876LI',
                       'G18923CQ','G23487NW','G26527XD','G31649CJ','G34368WA','G34723AD','G38409FH',
                       'G40748TG','G50188XO','G51941HK'):

            # WURCS hxh monosaccharide ambiguity...

            gctstr2 = glycoct_parser.toStr(g2)

            for l1,l2 in zip(gctstr.split('\n'),gctstr2.split('\n')):
                char = " "
                if l1 != l2:
                    char = "!"
                print("%-40s"%(l1[:40],),char,"%-40s"%(l2[:40],),file=sys.stderr)

            print(sum(1 for _ in g.all_nodes(subst=True)))
            print(sum(1 for _ in g2.all_nodes(subst=True)))
            for m1,m2 in zip(g.all_nodes(subst=True),g2.all_nodes(subst=True)):
                print(str(m1))
                print(str(m2))

            assert g.equals(g2)

        if g3:
            print(glycoct_parser.toStr(g3),file=sys.stderr)
            # assert g.equals(g3)

        j += 1

        if typ == "Saccharide":

            print(level(g),file=sys.stderr)
            
            assert level(g) == typ

            if gtopo == None or gcomp == None:
                continue

            # print "---------- sacc:topo ----------"

            # print glycoct_parser.toStr(gtopo)
            # print glycoct_parser.toStr(topology(g))

            print(glycoct_parser.toStr(gtopo),file=sys.stderr)
            print(glycoct_parser.toStr(topology(g)),file=sys.stderr)

            showdiff(gtopo,topology(g))
            
            assert g.equals(g1)
            assert not g.equals(gtopo)
            assert gtopo.equals(topology(g))

            print("---------- sacc:comp ----------")

            print(glycoct_parser.toStr(gcomp),file=sys.stderr)
            print(glycoct_parser.toStr(composition(g)),file=sys.stderr)

            assert not g.equals(gcomp)
            assert gcomp.equals(composition(g))
            assert g.equals(g1)

        elif typ == "Topology":

            assert level(g) == typ

            if gcomp == None:
                continue

            assert g.equals(topology(g))
            assert not g.equals(gcomp)
            assert gcomp.equals(composition(g))
            assert g.equals(g1)

        elif typ == "Composition":

            assert g.equals(composition(g))
            assert g.equals(g1)
