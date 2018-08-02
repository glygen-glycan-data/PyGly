
from Monosaccharide import Anomer, Stem, SuperClass, Config

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
            for l in m.links(instantiated_only=False):
                if l.parent_pos() == set([1]) and l.child_pos() == set([1]) and m.anomer() not in (Anomer.alpha, Anomer.beta):
                    l.set_parent_pos(None)
                if l.parent_pos() != set([1]):
                    l.set_parent_pos(None)
                if l.child_pos() != set([1]) and l.child().superclass() in (SuperClass.HEX,):
                    l.set_child_pos(None)
                if l.undetermined() and l.child_pos() == set([1]) and l.child().superclass() in (SuperClass.NON,SuperClass.OCT):
                    l.set_child_pos(None)
            if m.anomer() != Anomer.uncyclized:
                m.set_anomer(Anomer.missing)
        undets = list(g.undetermined_roots())
        g.set_undetermined(undets)

class Composition(Manipulation):

    topo = Topology()

    def manip(self,g):
        self.topo.manip(g)
        undets = []
        for m in g.all_nodes():
            undets.append(m)
        for m in undets:
            m.set_connected(False)
            m.clear_links()
            m.clear_parent_links()
            if "%s:%s"%(m.ring_start(),m.ring_end()) == "0:0":
                continue
            if m.anomer() == Anomer.uncyclized:
                continue
            m.set_ring_start(None)
            m.set_ring_end(None)
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

if __name__ == "__main__":

    import sys, csv, zipfile, traceback

    from GlycanFormatter import WURCS20Format, GlycoCTFormat, WURCS20ParseError, UnsupportedLinkError,  \
				CircularError, GlycoCTParseError, MonoOrderLinkError
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

        print >>sys.stderr, "\t".join(map(str,[i,j,acc,typ,topo,comp,bcomp]))

        try:
            gtopo = None; gcomp = None; 
            g = wurcs_parser.toGlycan(zf.read(acc+'.txt'))
            if topo:
                gtopo = wurcs_parser.toGlycan(zf.read(topo+'.txt'))
            if comp:
                gcomp = wurcs_parser.toGlycan(zf.read(comp+'.txt'))
        except KeyError:
	    print >>sys.stderr, "Skip: Accession not in zip file"
            # accession not in zip file
            continue
        except UnsupportedLinkError:
	    print >>sys.stderr, "Skip: Unsupported link error"
            # unexpected link structure
            continue
	except MonoOrderLinkError:
	    print >>sys.stderr, "Skip: Unsupported monosaccharide order"
	    # link order issue
	    continue
        except CircularError:
	    print >>sys.stderr, "Skip: Unsupported circular glycan"
            # dodge this case...
            continue
        except UnsupportedMonoError, e:
            print >>sys.stderr, "Skip:",e.message
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
	    print >>sys.stderr, "Skip: (1+1) two monosaccharide issue"
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

        print >>sys.stderr, glycoct_parser.toStr(g)

        g1 = g.clone()

        assert g.equals(g1)

	gctstr = glycoct_parser.toStr(g)

	g2 = glycoct_parser.toGlycan(gctstr)

	if acc not in ('G00727XI','G08545ST','G09485LV','G11795CV','G13902SG','G14164XT','G17876LI',
                       'G18923CQ','G23487NW','G26527XD','G31649CJ','G34368WA','G34723AD','G38409FH',
                       'G40748TG','G50188XO','G51941HK'):

	    # WURCS hxh monosaccharide ambiguity...

            gctstr2 = glycoct_parser.toStr(g2)

            for l1,l2 in zip(gctstr.split('\n'),gctstr2.split('\n')):
                char = " "
                if l1 != l2:
                    char = "!"
                print >>sys.stderr, "%-40s"%(l1[:40],),char,"%-40s"%(l2[:40],)

	    assert g.equals(g2)

	if g3:
            print >>sys.stderr, glycoct_parser.toStr(g3)
	    # assert g.equals(g3)

        j += 1

        if typ == "Saccharide":

	    print >>sys.stderr, level(g)
	    
            assert level(g) == typ

            if gtopo == None or gcomp == None:
                continue

            # print "---------- sacc:topo ----------"

            # print glycoct_parser.toStr(gtopo)
            # print glycoct_parser.toStr(topology(g))

            print >>sys.stderr, glycoct_parser.toStr(gtopo)
            print >>sys.stderr, glycoct_parser.toStr(topology(g))

            assert g.equals(g1)
            assert not g.equals(gtopo)
            assert gtopo.equals(topology(g))

            # print "---------- sacc:comp ----------"

            print >>sys.stderr, glycoct_parser.toStr(gcomp)
            print >>sys.stderr, glycoct_parser.toStr(composition(g))

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
