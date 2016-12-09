# -*- coding: utf-8 -*-
import enum
import itertools

import networkx

from camsa.core.data_structures import ScaffoldAssemblyGraph


class Conflicts(enum.Enum):
    non_conflicted = 0
    semi_conflicted = 1
    conflicted = 2


def edges_conflict(edge1, edge2):
    u1, v1 = edge1
    u2, v2 = edge2
    return (u1 in [u2, v2] and v1 not in [u2, v2]) or (u1 not in [u2, v2] and v1 in [u2, v2])


def get_conflict_type(ap1, ap2):
    conflicted_cnt = 0
    non_conflicted_cnt = 0
    for edge1, edge2 in itertools.product(ap1.get_edges(), ap2.get_edges()):
        if edges_conflict(edge1=edge1, edge2=edge2):
            conflicted_cnt += 1
        else:
            non_conflicted_cnt += 1
    if conflicted_cnt == 0:
        return Conflicts.non_conflicted
    elif conflicted_cnt > 0 and non_conflicted_cnt > 0:
        return Conflicts.semi_conflicted
    else:
        return Conflicts.conflicted


def get_conflicting_edges(sag, ap_edge):
    result = []
    for vertex in ap_edge:
        conflicting_edges = [edge for edge in sag.graph.edges(nbunch=vertex, data=True)
                             if edges_conflict(edge1=(edge[0], edge[1]), edge2=ap_edge)]
        result.extend(conflicting_edges)
    return result


def get_conflicting_assembly_points(sag, assembly_point, assembly_points_by_ids):
    conflicting_assembly_points_ids = []
    for edge in assembly_point.get_edges():
        conflicting_edges = get_conflicting_edges(sag=sag, ap_edge=edge)
        ids = [data["ap_id"] for (u, v, data) in conflicting_edges]
        conflicting_assembly_points_ids.extend(ids)
    return [assembly_points_by_ids[ap_id] for ap_id in set(conflicting_assembly_points_ids) if ap_id != assembly_point.self_id]


def update_conflicts(ap1, ap2, in_ap1, out_ap1, in_ap2, out_ap2):
    in_sources = set(ap1.sources) & set(ap2.sources)
    ap1_out_sources = set(ap2.sources) - set(ap1.sources)
    ap2_out_sources = set(ap1.sources) - set(ap2.sources)
    for source in in_sources:
        in_ap1[source].add(ap2.self_id)
        in_ap2[source].add(ap1.self_id)
    for source in ap1_out_sources:
        out_ap1[source].add(ap2.self_id)
    for source in ap2_out_sources:
        out_ap2[source].add(ap1.self_id)


def update_assembly_points_as_conflicted(ap1, ap2):
    in_ap1_destination = ap1.in_conflicted
    out_ap1_destination = ap1.out_conflicted
    in_ap2_destination = ap2.in_conflicted
    out_ap2_destination = ap2.out_conflicted
    update_conflicts(ap1=ap1, ap2=ap2,
                     in_ap1=in_ap1_destination, out_ap1=out_ap1_destination,
                     in_ap2=in_ap2_destination, out_ap2=out_ap2_destination)


def update_assembly_points_as_semi_conflicted(ap1, ap2):
    in_ap1_destination = ap1.in_semi_conflicted
    out_ap1_destination = ap1.out_semi_conflicted
    in_ap2_destination = ap2.in_semi_conflicted
    out_ap2_destination = ap2.out_semi_conflicted
    update_conflicts(ap1=ap1, ap2=ap2,
                     in_ap1=in_ap1_destination, out_ap1=out_ap1_destination,
                     in_ap2=in_ap2_destination, out_ap2=out_ap2_destination)


def compute_and_update_assembly_points_conflicts(assembly_points_by_ids):
    """

    :param assembly_points_by_ids: a dictionary, where key is the assembly point id, and value is merged assembly points
    :return:
    """
    sag = ScaffoldAssemblyGraph.from_assembly_points(assembly_points=assembly_points_by_ids.values(),
                                                     add_scaffold_edges=False)
    for ap in assembly_points_by_ids.values():
        conflicted_assembly_points = get_conflicting_assembly_points(sag=sag, assembly_point=ap,
                                                                     assembly_points_by_ids=assembly_points_by_ids)
        for c_ap in conflicted_assembly_points:
            conflict_type = get_conflict_type(ap1=ap, ap2=c_ap)
            assert conflict_type != Conflicts.non_conflicted
            if conflict_type == Conflicts.conflicted:
                update_assembly_points_as_conflicted(ap1=ap, ap2=c_ap)
            elif conflict_type == Conflicts.semi_conflicted:
                update_assembly_points_as_semi_conflicted(ap1=ap, ap2=c_ap)
