import itertools
from pygly.Monosaccharide import Linkage

class motif_search_dev():

    def all_nodes(self,glycan_root):
        nodes = [glycan_root]
        self.append_mono(glycan_root,nodes)

    def append_mono(self,mono, nodeslist):
        links = list(mono.links())
        if len(links) != 0:
            for link in links:
                child = link.child()
                nodeslist.append(child)
                self.append_mono(child, nodeslist)

    def mono_compatible_one_way(self, mg, mm):
        if (mm._anomer) and mg._anomer != mm._anomer:
            return False
        if mg._config != mm._config:
            return False
        if mg._stem != mm._stem:
            return False
        if mg._superclass != mm._superclass:
            return False
        if mg._ring_start != mm._ring_start:
            return False
        if mg._ring_end != mm._ring_end:
            return False
        if mg._mods != mm._mods:
            return False

        subsg = mg._substituent_links
        subsm = mm._substituent_links
        if len(subsg) != len(subsm):
            return False
        for newsubsm in itertools.permutations(subsm):
            flag_djkswjw = True
            for subg, subm in zip(subsg, newsubsm):
                flag_djkswjw = self.link_compatible_one_way(subg, subm) and subg.child().equals(subm.child())
                if not flag_djkswjw:
                    break
            if flag_djkswjw:
                break
        if not flag_djkswjw:
            return False
        return True

    def mono_compatible_one_way2(self, mg, mm):
        if (mm._anomer) and mg._anomer != mm._anomer:
            return False
        if (mm._config) and mg._config != mm._config:
            return False
        if (mm._stem) and mg._stem != mm._stem:
            return False
        if mg._superclass != mm._superclass:
            return False
        if (mm._ring_start) and mg._ring_start != mm._ring_start:
            return False
        if (mm._ring_end) and mg._ring_end != mm._ring_end:
            return False
        if mg._mods != mm._mods:
            return False

        subsg = mg._substituent_links
        subsm = mm._substituent_links
        if len(subsg) != len(subsm):
            return False
        for newsubsm in itertools.permutations(subsm):
            flag_djkswjw = True
            for subg, subm in zip(subsg, newsubsm):
                flag_djkswjw = self.link_compatible_one_way(subg, subm) and subg.child().equals(subm.child())
                if not flag_djkswjw:
                    break
            if flag_djkswjw:
                break
        if not flag_djkswjw:
            return False
        return True

    def mono_equals(self, mg, mm):
        return mg.equals(mm)

    def link_compatible_one_way(self, lg, lm):
        if lg._parent_type != lm._parent_type:
            return False
        if lg._child_type != lm._child_type:
            return False
        if lg._child_pos != lm._child_pos:
            return False
        ppg = Linkage.parentpos2set(lg)
        ppm = Linkage.parentpos2set(lm)
        if (ppg.issubset(ppm) and len(ppg) != 0) or len(ppm) == 0:
            return True
        return False

    def link_compatible_one_way2(self, lg, lm):
        if lg._parent_type != lm._parent_type:
            return False
        if lg._child_type != lm._child_type:
            return False
        if lg._child_pos != lm._child_pos:
            return False
        ppg = Linkage.parentpos2set(lg)
        ppm = Linkage.parentpos2set(lm)
        if len(ppm) == 0 or len(ppg) == 0 or len(ppm.union(ppg)) < (len(ppg) + len(ppm)):
            return True
        return False

    def link_equals(self, lg, lm):
        return lg.equals(lm)

    def mono_check(self,a,b):
        raise NotImplemented()

    def search(self, any_node, motif_root_sub):
        linksg = list(any_node.links())
        linksm = list(motif_root_sub.links())
        if self.mono_check(any_node,motif_root_sub) and (len(linksm) <= len(linksg)):
            if len(linksm) == 0:
                return True
            else:
                link1 = linksm[:]
                for link2 in itertools.permutations(linksg):
                    for l1, l2 in zip(link1, link2):
                        equallinks = self.link_check(l2,l1)
                        equalmono = self.search(l2.child(), l1.child())
                        flag = equallinks and equalmono
                        if not flag:
                            break
                    if flag:
                        return True
                return False
        else:
            return False

    def match_reducing_end(self, glycan_root, motif_root):
        return self.search(glycan_root,motif_root)

    def match_else_where(self, glycan_root, motif_root):
        nodes_except_for_root = []
        self.append_mono(glycan_root, nodes_except_for_root)
        for node in nodes_except_for_root:
            if self.search(node, motif_root):
                return True
        return False

class motif_search_accept_wildcard(motif_search_dev):

    def mono_check(self, mg, mm):
        return self.mono_compatible_one_way2(mg, mm)

    def link_check(self, lg, lm):
        return self.link_compatible_one_way2(lg, lm)

    def match(self,g,m):
        return self.match_reducing_end(g,m) or self.match_else_where(g,m)

class motif_search_strict(motif_search_dev):

    def mono_check(self, mg, mm):
        return self.mono_equals(mg, mm)

    def link_check(self, lg, lm):
        return self.link_equals(lg, lm)

    def match(self,g,m):
        return self.match_reducing_end(g, m) or self.match_else_where(g, m)

if __name__ == "__main__":
    pass