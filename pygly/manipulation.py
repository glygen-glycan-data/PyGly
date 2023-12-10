
from __future__ import print_function

from . Monosaccharide import Anomer, Stem, SuperClass, Config, Substituent

class Manipulation(object):

    def __call__(self,g):
        return self.make(g)

    def make(self,g):
        g1 = g.clone()
        self.manip(g1)
        return g1

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
