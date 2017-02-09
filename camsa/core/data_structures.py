# -*- coding: utf-8 -*-
import itertools
import json
from collections import defaultdict

import networkx
import six
import sys


class AssemblyPoint(object):
    def __init__(self, seq1, seq2, seq1_or, seq2_or, sources, cw=None, parent_id=None, children_ids=None, self_id=None,
                 gap_size=None):
        self.seq1 = seq1
        self.seq2 = seq2
        self.seq1_or = seq1_or
        self.seq2_or = seq2_or
        self.sources = sorted(sources)
        self.cw = cw
        self.gap_size = gap_size
        self.self_id = self_id
        self.parent_id = parent_id
        self.children_ids = children_ids if children_ids is not None else []
        #### comparative metrics
        self.participates_in_merged = False
        self.in_conflicted = defaultdict(set)
        self.in_semi_conflicted = defaultdict(set)
        self.out_conflicted = defaultdict(set)
        self.out_semi_conflicted = defaultdict(set)
        self.seq1_par_or = None
        self.seq2_par_or = None
        #### ref metrics
        self.ref_metrics = RefMetrics()

    @property
    def orientation_as_word(self):
        if self.is_unoriented:
            return "U"
        if self.is_semi_oriented:
            return "SO"
        return "O"

    @property
    def is_unoriented(self):
        return [self.seq1_or, self.seq2_or].count("?") == 2

    @property
    def is_semi_oriented(self):
        return [self.seq1_or, self.seq2_or].count("?") == 1

    @property
    def is_oriented(self):
        return [self.seq1_or, self.seq2_or].count("?") == 0

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

    def get_edges(self, weight=False, sort=True):
        result = []
        for ap in self.get_all_all_possible_realizations():
            head = ap.seq1 + ("h" if ap.seq1_or == "+" else "t")
            tail = ap.seq2 + ("t" if ap.seq2_or == "+" else "h")
            hv = head
            tv = tail
            if sort:
                hv, tv = (hv, tv) if hv < tv else (tv, hv)
            if weight:
                result.append((hv, tv, self.cw))
            else:
                result.append((hv, tv))
        return result

    def get_all_all_possible_realizations(self):
        result = []
        if self.seq1_or == "?":
            seq1_choices = ["+", "-"]
        else:
            seq1_choices = [self.seq1_or]
        if self.seq2_or == "?":
            seq2_choices = ["+", "-"]
        else:
            seq2_choices = [self.seq2_or]
        for ctg1_or, ctg2_or in itertools.product(seq1_choices, seq2_choices):
            result.append(AssemblyPoint(seq1=self.seq1, seq2=self.seq2, seq1_or=ctg1_or, seq2_or=ctg2_or,
                                        sources=self.sources))
        return result

    @property
    def is_ambiguous(self):
        return self.seq1_or == "?" or self.seq2_or == "?"


class APFieldOutExtractorConverter(object):
    def __init__(self, field_name, converter_name):
        self.field_name = field_name
        self.converter = AP_FIELD_CONVERTERS[converter_name]

    def extract_field_value_str(self, ap):
        return self.converter.convert(getattr(ap, self.field_name))


class APFieldConverter(object):
    @staticmethod
    def convert(field_value, **kwargs):
        raise NotImplemented("Converter must be field/area specific")


class APFieldConverterStr(APFieldConverter):
    @staticmethod
    def convert(field_value, **kwargs):
        return str(field_value)


class APFieldConverterQuestionIfNone(APFieldConverterStr):
    @staticmethod
    def convert(field_value, **kwargs):
        if field_value is None:
            field_value = "?"
        return super(APFieldConverterQuestionIfNone, APFieldConverterQuestionIfNone).convert(field_value=field_value, **kwargs)


class APFieldConverterBooleanToInt(APFieldConverterQuestionIfNone):
    @staticmethod
    def convert(field_value, **kwargs):
        if field_value in [True, False]:
            field_value = int(field_value)
        return super(APFieldConverterBooleanToInt, APFieldConverterBooleanToInt).convert(field_value=field_value, **kwargs)


IdFieldConverter = APFieldConverterQuestionIfNone


class IterableFieldConverter(APFieldConverterQuestionIfNone):
    @staticmethod
    def convert(field_value, **kwargs):
        result = []
        treat_as_set = kwargs.get("as_set", True)
        for value in field_value:
            result.append(super(IterableFieldConverter, IterableFieldConverter).convert(field_value=value, **kwargs))
        if treat_as_set:
            result = set(result)
        result = sorted(result)
        treat_empty_as_zero = kwargs.get("empty_as_zero", True)
        if treat_empty_as_zero and len(result) == 0:
            return "0"
        separator = kwargs.get("intra_separator", ",")
        return separator.join(result)


class ConflictFieldConverter(APFieldConverter):
    @staticmethod
    def convert(field_value):
        conflicted_ids = [ap_id for conflict_assembly in field_value.values() for ap_id in conflict_assembly]
        return IterableFieldConverter.convert(field_value=conflicted_ids)


class RefMetrics(object):
    def __init__(self):
        #### correct
        self.present = False
        self.present_ref_ids = []
        #### incorrect
        # the worst error
        self.interchromosomal = False
        # the best error
        self.tranlocation = False
        self.translocation_path = []
        # can be combined with translocation, or can be on its own
        self.reversal = False
        # a subtype of reversal, better than reversal. possible for semi/un-oriented assembly points only
        self.semi_reversal = False
        # a best available reading of the reference genome for this particular assembly point
        self.best_ref_reading_id = None


class MergedScaffoldAssemblyGraph(object):
    def __init__(self):
        self.graph = networkx.Graph()

    def add_edge(self, u, v, weight):
        if self.graph.has_edge(u=u, v=v):
            self.graph[u][v]['weight'] += weight
        else:
            self.graph.add_edge(u=u, v=v, weight=weight)

    def get_maximal_non_conflicting_assembly_graph(self):
        max_matching = networkx.max_weight_matching(G=self.graph)
        seen = set()
        edges = []
        for key, value in max_matching.items():
            if key not in seen and value not in seen:
                edges.append((key, value))
                seen.add(key)
                seen.add(value)
        result = networkx.Graph()
        result.add_edges_from(edges)
        return result

    def edges(self, weight=True):
        if weight:
            return [(u, v, w) if u < v else (v, u, w) for (u, v, w) in self.graph.edges(data='weight')]
        else:
            return [(u, v) if u < v else (v, u) for (u, v) in self.graph.edges()]

    def remove_edges_with_low_cw(self, cw_threshold=0.0):
        low_weight_edges = [(u, v) for (u, v, d) in self.graph.edges(data=True) if d['weight'] < cw_threshold]
        self.graph.remove_edges_from(low_weight_edges)


class ScaffoldAssemblyGraph(object):
    def __init__(self):
        self.graph = networkx.MultiGraph()

    def add_edge(self, u, v, **kwargs):
        self.graph.add_edge(u=u, v=v, attr_dict=kwargs)

    @classmethod
    def from_assembly_points(cls, assembly_points, reference=False, add_scaffold_edges=True):
        result = cls()
        scaffold_edges = get_scaffold_edges(assembly_points=assembly_points)
        if add_scaffold_edges:
            for u, v in scaffold_edges:
                result.add_edge(u=u, v=v, edge_type="scaffold")
        for ap in assembly_points:
            for (u, v, weight) in ap.get_edges(sort=True, weight=True):
                edge_type = "ref_assembly" if reference else "assembly"
                result.add_edge(u=u, v=v, weight=weight, edge_type=edge_type, ap_id=ap.self_id, ap_type=ap.orientation_as_word)
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
        return len([ap for ap in self.aps if ap.is_non_conflicted])


def assign_ids_to_assembly_points(assembly_points, id_prefix="", id_generator=None):
    assembly_points_by_ids = {}
    if id_generator is None:
        id_generator = itertools.count()
    for assembly_point in assembly_points:
        assembly_point.self_id = id_prefix + str(six.next(id_generator))
        assembly_points_by_ids[assembly_point.self_id] = assembly_point
    return assembly_points_by_ids


def merge_assembly_points(assembly_points_by_source):
    unique_assembly_points = defaultdict(list)

    for assembly_points in assembly_points_by_source.values():
        for assembly_point in assembly_points:
            seq1, seq2 = assembly_point.seq1, assembly_point.seq2
            seq1_or, seq2_or = assembly_point.seq1_or, assembly_point.seq2_or
            if seq1 < seq2:
                entry = (seq1, seq1_or, seq2, seq2_or)
            else:
                entry = (seq2, inverse_orientation(seq2_or), seq1, inverse_orientation(seq1_or))
            unique_assembly_points[entry].append(assembly_point)
    result = []
    for (seq1, seq1_or, seq2, seq2_or), children in unique_assembly_points.items():
        sources = sorted(set(source for ap in children for source in ap.sources))
        weight = sum(ap.cw for ap in children)
        children_ids = [ap.self_id for ap in children]
        merged_assembly_point = AssemblyPoint(seq1=seq1, seq2=seq2, seq1_or=seq1_or, seq2_or=seq2_or,
                                              sources=sources, cw=weight, children_ids=children_ids)
        result.append(merged_assembly_point)
    return result


def assign_parents_to_children(children_assembly_points_by_ids, parent_assembly_points_by_ids):
    for p_assembly_point in parent_assembly_points_by_ids.values():
        for child_id in p_assembly_point.children_ids:
            children_assembly_points_by_ids[child_id].parent_id = p_assembly_point.self_id


def inverse_orientation(orientation):
    if orientation in ("+", "-"):
        return "+" if orientation == "-" else "-"
    return orientation


AP_FIELD_CONVERTERS = {
    "str_raw": APFieldConverterStr,
    "str": APFieldConverterQuestionIfNone,
    "bool": APFieldConverterBooleanToInt,
    "conflict": ConflictFieldConverter,
    "id": IdFieldConverter,
    "iter": IterableFieldConverter,
}


def get_scaffold_edges(assembly_points):
    unique_scaffolds = set()
    for ap in assembly_points:
        unique_scaffolds.add(ap.seq1)
        unique_scaffolds.add(ap.seq2)
    return [(name + "t", name + "h") for name in unique_scaffolds]


def to_json(value):
    if isinstance(value, Sequence):
        return json.dumps({"seq_id": value.name, "length": value.length})
    result = {}
    for key, values in value.items():
        if isinstance(values, str):
            result[key] = values
        else:
            result[key] = sorted(values)
    return json.dumps(result)


class Sequence(object):
    def __init__(self, name, length=None, parent_seq_id=None, start=None, end=None, strand="+", annotation=None, seq_group_id="Default"):
        self.name = name
        self._length = length
        self._parent_seq_id = parent_seq_id
        self._start = start
        self._end = end
        self._seq_group_id = seq_group_id
        self.strand = strand
        self._annotation = annotation

    @property
    def length(self):
        if self._length is not None:
            return self._length
        if self.start is not None and self.end is not None:
            return self.end - self.start
        return None

    @property
    def parent_seq_id(self):
        if self._parent_seq_id is None:
            return self.name
        return self._parent_seq_id

    @property
    def start(self):
        if self._start is None:
            return 0
        return int(self._start)

    @property
    def end(self):
        if self._end is None:
            return sys.maxsize
        return int(self._end)

    @property
    def seq_group_id(self):
        return str(self._seq_group_id)

    @property
    def annotation(self):
        return str(self._annotation)
