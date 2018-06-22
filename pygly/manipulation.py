
from Monosaccharide import Anomer, Stem

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
            if m.anomer() != Anomer.uncyclized:
                m.set_anomer(Anomer.missing)
            for l in m.links(instantiated_only=False):
                if l.parentpos2set() == set([1]):
                    continue
                l.set_parent_pos(None)
                l.set_parent_pos2(None)
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

class LevelSniffer(object):
    topology = Topology()
    composition = Composition()
    basecomposition = BaseComposition()
    
    def __call__(self,g):
        return self.sniff(g)

    def sniff(self,g):
        if g.hasroot():
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

    import sys, csv, zipfile

    from GlycanFormatter import WURCS20Format, GlycoCTFormat, WURCS20ParseError

    wurcs_parser = WURCS20Format()
    glycoct_parser = GlycoCTFormat()

    topology = Topology()
    composition = Composition()
    basecomposition = BaseComposition()
    level = LevelSniffer()
    
    zf = zipfile.ZipFile(sys.argv[1])
    
    reader = csv.DictReader(sys.stdin,dialect='excel-tab')

    j = 0
    for i,d in enumerate(reader):
        acc = d['GlyTouCanAccession']
        typ = d['GlyTouCanType']
        topo = d.get('Topology')
        comp = d.get('Composition')
        bcomp = d.get('BaseComposition')

        print "\t".join(map(str,[i,j,acc,typ,topo,comp,bcomp]))

        try:
            gtopo = None; gcomp = None; 
            g = wurcs_parser.toGlycan(zf.read(acc+'.txt'))
            if topo:
                gtopo = wurcs_parser.toGlycan(zf.read(topo+'.txt'))
            if comp:
                gcomp = wurcs_parser.toGlycan(zf.read(comp+'.txt'))
        except (WURCS20ParseError,KeyError):
            continue

        print glycoct_parser.toStr(g)

        g1 = g.clone()

        assert g.equals(g1)
            
        j += 1

        if typ == "Saccharide":

	    print level(g)
            print glycoct_parser.toStr(topology(g))
	    
            assert level(g) == typ

            if gtopo == None or gcomp == None:
                continue

            # print "---------- sacc:topo ----------"

            # print glycoct_parser.toStr(gtopo)
            # print glycoct_parser.toStr(topology(g))

            assert g.equals(g1)
            assert not g.equals(gtopo)
            assert gtopo.equals(topology(g))

            # print "---------- sacc:comp ----------"

            # print glycoct_parser.toStr(gcomp)
            # print glycoct_parser.toStr(composition(g))

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

            assert level(g) == typ

            assert g.equals(composition(g))
            assert g.equals(g1)
