# -*- coding: utf-8 -*-
import itertools

import networkx
from bg.edge import BGEdge
from bg.genome import BGGenome
from bg.multicolor import Multicolor
from bg.vertices import TaggedBlockVertex


class AssemblyPoint(object):
    def __init__(self, ctg1, ctg2, ctg1_or, ctg2_or, sources, cw=None, parent_id=None, children_id=None, self_id=None):
        self.contig_1 = ctg1
        self.contig_2 = ctg2
        self.contig_1_orientation = ctg1_or
        self.contig_2_orientation = ctg2_or
        self.sources = sorted(sources)
        self.participates_in_max_non_conflicting_assembly = False
        self.in_conflicted = []
        self.in_semi_conflicted = []
        self.out_conflicted = []
        self.out_semi_conflicted = []
        self.participation_ctg1_or = None
        self.participation_ctg2_or = None
        self.cw = cw
        self.self_id = self_id
        self.parent_id = parent_id
        self.children_id = children_id if children_id is not None else []

    @property
    def orientation_as_word(self):
        if self.is_unoriented: return "U"
        if self.is_semi_oriented: return "SO"
        return "O"

    def get_bg_edges(self):
        result = []
        for ap in self.get_all_all_possible_representations():
            head = ap.contig_1 + ("h" if ap.contig_1_orientation == "+" else "t")
            tail = ap.contig_2 + ("t" if ap.contig_2_orientation == "+" else "h")
            hv = TaggedBlockVertex(name=head)
            tv = TaggedBlockVertex(name=tail)
            result.append(
                BGEdge(vertex1=hv, vertex2=tv, multicolor=Multicolor(*[BGGenome(name) for name in self.sources])))
        return result

    @property
    def is_unoriented(self):
        return [self.contig_1_orientation, self.contig_2_orientation].count("?") == 2

    @property
    def is_semi_oriented(self):
        return [self.contig_1_orientation, self.contig_2_orientation].count("?") == 1

    @property
    def is_oriented(self):
        return [self.contig_1_orientation, self.contig_2_orientation].count("?") == 0

    def is_in_semi_conflicted_for(self, source_name):
        return source_name in self.in_semi_conflicted

    def is_in_conflicted_for(self, source_name):
        return source_name in self.in_conflicted

    def is_out_semi_conflicted_for(self, source_name):
        return len(self.out_semi_conflicted) > 0

    def is_out_conflicted_for(self, source_name):
        return len(self.out_conflicted) > 0

    @property
    def is_out_conflicted(self):
        return len(self.out_conflicted) > 0

    @property
    def is_out_semi_conflicted(self):
        return len(self.out_semi_conflicted) > 0

    @property
    def is_non_conflicted(self):
        return len(self.in_conflicted) == 0 and len(self.out_conflicted) == 0 and len(self.in_semi_conflicted) == 0 and len(
            self.out_semi_conflicted) == 0

    def get_edges(self):
        result = []
        for ap in self.get_all_all_possible_representations():
            head = ap.contig_1 + ("h" if ap.contig_1_orientation == "+" else "t")
            tail = ap.contig_2 + ("t" if ap.contig_2_orientation == "+" else "h")
            hv = TaggedBlockVertex(name=head)
            tv = TaggedBlockVertex(name=tail)
            # hv = head
            # tv = tail
            result.append((hv, tv))
        return result

    def get_all_all_possible_representations(self):
        result = []
        if self.contig_1_orientation == "?":
            ctg_1_choices = ["+", "-"]
        else:
            ctg_1_choices = [self.contig_1_orientation]
        if self.contig_2_orientation == "?":
            ctg_2_choices = ["+", "-"]
        else:
            ctg_2_choices = [self.contig_2_orientation]
        for ctg1_or, ctg2_or in itertools.product(ctg_1_choices, ctg_2_choices):
            result.append(AssemblyPoint(ctg1=self.contig_1, ctg2=self.contig_2, ctg1_or=ctg1_or, ctg2_or=ctg2_or,
                                        sources=self.sources))
        return result

    @property
    def is_ambiguous(self):
        return self.contig_1_orientation == "?" or self.contig_2_orientation == "?"


class AssemblyGraph(object):
    def __init__(self):
        self.graph = networkx.Graph()

    def add_edge(self, u, v, weight=1):
        if self.graph.has_edge(u=u, v=v):
            self.graph[u][v]['weight'] += weight
        else:
            self.graph.add_edge(u=u, v=v, weight=weight)

    def get_maximal_non_conflicting_assembly_graph(self):
        max_matching = networkx.max_weight_matching(G=self.graph)
        seen = set()
        edges = []
        sum = 0
        for key, value in max_matching.items():
            if key not in seen and value not in seen:
                edges.append((key, value))
                sum += self.graph[key][value]["weight"]
                if 0.57 < self.graph[key][value]["weight"] < 0.59:
                    print(key, value, self.graph[key][value]["weight"])
                seen.add(key)
                seen.add(value)
        result = networkx.Graph()
        print(sum)
        result.add_edges_from(edges)
        return result


class Assembly(object):
    def __init__(self, name, aps):
        self.name = name
        self.aps = aps
        self.total_cnt = 0
        self.max_non_conflicting_aps = []

    @property
    def entries_names(self):
        return self.get_all_entries_names()

    def get_all_entries_names(self):
        return sorted(set(name for ap in self.aps for name in ap.sources))

    def sort_aps(self):
        self.aps = sorted(self.aps,
                          key=lambda ap: (ap.contig_1, ap.contig_1_orientation, ap.contig_2, ap.contig_2_orientation))

    @property
    def unoriented_ap_cnt(self):
        return sum(ap.is_unoriented for ap in self.aps)

    @property
    def semi_oriented_aps_cnt(self):
        return sum(ap.is_semi_oriented for ap in self.aps)

    @property
    def oriented_aps_cnt(self):
        return sum(ap.is_oriented for ap in self.aps)

    @property
    def in_conflicted_cnt(self):
        cnt = 0
        for ap in self.aps:
            if isinstance(self.name, tuple):
                to_check = self.name
            else:
                to_check = [self.name]
            for name in to_check:
                if ap.is_in_conflicted_for(name):
                    cnt += 1
                    break
        return cnt

    @property
    def in_semi_conflicted_cnt(self):
        cnt = 0
        for ap in self.aps:
            if isinstance(self.name, tuple):
                to_check = self.name
            else:
                to_check = [self.name]
            for name in to_check:
                if ap.is_in_semi_conflicted_for(name):
                    cnt += 1
                    break
        return cnt

    @property
    def out_conflicted_cnt(self):
        return sum(ap.is_out_conflicted for ap in self.aps)

    @property
    def out_semi_conflicted_cnt(self):
        return sum(ap.is_out_semi_conflicted for ap in self.aps)

    @property
    def non_conflicted_cnt(self):
        return len([ap for ap in self.aps if (not ap.is_in_semi_conflicted_for(self.name)) and
                    (not ap.is_in_conflicted_for(self.name)) and
                    (not ap.is_out_conflicted) and
                    (not ap.is_out_semi_conflicted)])
