#! /usr/bin/env python
# -*- coding: utf-8 -*-
import datetime
import os
import shutil
from argparse import ArgumentParser
from collections import defaultdict

import itertools

from data_structures import AssemblyPoint, AssemblyGraph, Assembly
from jinja2 import Template
import camsa_io as camsa_io


VERSION = "1.0.0"


def inverse_orientation(orientation):
    if orientation in ("+", "-"):
        return "+" if orientation == "-" else "-"
    return orientation


def delete_irregular_edges(breakpoint_graph):
    irregular_edges = [edge for edge in breakpoint_graph.edges() if edge.is_irregular_edge]
    for edge in irregular_edges:
        breakpoint_graph.delete_bgedge(bgedge=edge)


def delete_irregular_vertices(breakpoint_graph):
    irregular_vertices = [vertex for vertex in breakpoint_graph.nodes() if vertex.is_irregular_vertex]
    for vertex in irregular_vertices:
        breakpoint_graph.bg.remove_node(vertex)


def get_non_conflicting_vertices(breakpoint_graph, vertices):
    result = []
    for vertex in vertices:
        if len(list(breakpoint_graph.get_edges_by_vertex(vertex=vertex))) == 1:
            result.append(vertex)
    return result


def get_conflicting_vertices(breakpoint_graph, vertices):
    result = []
    for vertex in vertices:
        if len(list(breakpoint_graph.get_edges_by_vertex(vertex=vertex))) > 1:
            result.append(vertex)
    return result


def get_conflicted_assemblies(breakpoint_graph, conflict_vertices):
    result = Assembly(name=None, aps=[])
    for vertex in conflict_vertices:
        v1 = vertex
        for edge in breakpoint_graph.get_edges_by_vertex(vertex=vertex):
            v2 = edge.vertex1 if v1 == edge.vertex2 else edge.vertex2
            ctg1 = v1.block_name
            ctg1_or = "-" if v1.is_tail_vertex else "+"
            ctg2 = v2.block_name
            ctg2_or = "-" if v2.is_head_vertex else "+"
            sources = [genome.name for genome in edge.multicolor.colors]
            ap = AssemblyPoint(ctg1=ctg1, ctg2=ctg2, ctg1_or=ctg1_or, ctg2_or=ctg2_or, sources=sources)
            result.aps.append(ap)
            result.total_cnt += 1
    return result


def read_aps_as_pairs(stream, separator="\t", destination=None):
    if destination is None:
        result = defaultdict(list)
    else:
        result = destination
    for line in stream:
        data = line.split(separator)
        result[data[0]].append(
            AssemblyPoint(ctg1=data[1], ctg2=data[2], ctg1_or=data[3], ctg2_or=data[4], sources=[data[0]]))
    return result


def get_non_conflicting_assemblies(breakpoint_graph, non_conflicting_vertices):
    result = {}
    for vertex in non_conflicting_vertices:
        v1 = vertex
        edge = next(breakpoint_graph.get_edges_by_vertex(vertex=vertex))
        v2 = edge.vertex1 if v1 == edge.vertex2 else edge.vertex2
        ctg1 = v1.block_name
        ctg1_or = "-" if v1.is_tail_vertex else "+"
        ctg2 = v2.block_name
        ctg2_or = "-" if v2.is_head_vertex else "+"
        sources = [genome.name for genome in edge.multicolor.colors]
        if len(sources) == 1:
            key = sources[0]
        else:
            key = tuple(sorted(sources))
        if key not in result:
            result[key] = Assembly(name=key, aps=[])
        if isinstance(key, tuple):
            for value in key:
                if value not in result:
                    result[value] = Assembly(name=value, aps=[])
        result[key].aps.append(AssemblyPoint(ctg1=ctg1, ctg2=ctg2, ctg1_or=ctg1_or, ctg2_or=ctg2_or, sources=sources))
        result[key].total_cnt += 1
        if isinstance(key, tuple):
            for value in key:
                result[value].total_cnt += 1
    assemblies = [value for key, value in result.items() if isinstance(key, str)]
    assemblies_intersections = [value for key, value in result.items() if isinstance(key, tuple)]
    return assemblies, assemblies_intersections


def get_assembly_points(breakpoint_graph):
    result = defaultdict(list)
    for edge in breakpoint_graph.edges():
        organism = edge.multicolor.colors.pop().name
        ctg1 = edge.vertex1.block_name
        ctg2 = edge.vertex2.block_name
        ctg1_orientation = "+" if edge.vertex1.is_head_vertex else "-"
        ctg2_orientaiton = "+" if edge.vertex2.is_tail_vertex else "-"
        ap = AssemblyPoint(ctg1=ctg1, ctg2=ctg2, ctg1_or=ctg1_orientation, ctg2_or=ctg2_orientaiton, sources=[organism])
        result[organism].append(ap)
    return result


def merge_aps(aps, or_ap_by_id=None, m_ap_by_id=None):
    sources = defaultdict(list)
    cw = defaultdict(int)
    children = defaultdict(list)
    or_ap_by_id = or_ap_by_id if or_ap_by_id is not None else {}
    m_ap_by_id = m_ap_by_id if m_ap_by_id is not None else {}
    for organism, o_aps in aps.items():
        for cnt, ap in enumerate(o_aps):
            ap.self_id = "or_{id}".format(id=cnt)
            or_ap_by_id[ap.self_id] = ap
            ctg1, ctg2 = ap.contig_1, ap.contig_2
            ctg1_or, ctg2_or = ap.contig_1_orientation, ap.contig_2_orientation
            ctg1_or, ctg2_or = (ctg1_or, ctg2_or) if ctg1 < ctg2 else (
                inverse_orientation(ctg2_or), inverse_orientation(ctg1_or))
            ctg1, ctg2 = (ctg1, ctg2) if ctg1 < ctg2 else (ctg2, ctg1)
            entry = (ctg1, ctg2, ctg1_or, ctg2_or)
            cw[entry] += ap.cw
            sources[entry].extend(ap.sources)
            children[entry].append(ap.self_id)
    result = []
    for cnt, (entry, sources) in enumerate(sources.items()):
        children_ids = children[entry]
        ap = AssemblyPoint(ctg1=entry[0], ctg2=entry[1], ctg1_or=entry[2], ctg2_or=entry[3], sources=sources, cw=cw[entry],
                           children_id=children_ids, self_id="m_{id}".format(id=cnt))
        m_ap_by_id[ap.self_id] = ap
        for ap_id in ap.children_id:
            or_ap_by_id[ap_id].parent_id = ap.self_id
        result.append(ap)
    return result


def get_all_conflicted_aps(ap, all_aps):
    conflicts = []
    for u, v in ap.get_edges():
        c_aps = {c_ap for c_ap in all_aps[u] if c_ap != ap}
        conflicts.append(c_aps.union({c_ap for c_ap in all_aps[v] if c_ap != ap}))
    conflicted = set.intersection(*conflicts)
    semi_conflicted = {ap for c_aps_set in conflicts for ap in c_aps_set if ap not in conflicted}
    return conflicted, semi_conflicted


if __name__ == "__main__":
    full_description = """
        Sergey Aganezov & Max A. Alekseyev (c)
        Computational Biology Institute, The George Washington University.

        CAMSA is a tool for Comparative Analysis and Merging of Scaffold Assemblies.
        For more information relate to README.md file in the root of CAMSA distribution.

        With any questions, please, contact Sergey Aganezov [aganezov(at)gwu.edu].
        """
    parser = ArgumentParser(description=full_description)
    parser.add_argument("input", nargs="+")
    parser.add_argument("-f", "--output-file", dest="output_report_file", default="report.html")
    parser.add_argument("--default-exact-cw", dest="cw_exact", type=float, default=1.0)
    parser.add_argument("--default-prob-cw", dest="cw_prob", type=float, default=0.9)
    parser.add_argument("--minimum-cw-threshold", dest="min_cw", type=float, default=0.0)
    parser.add_argument("--version", action="version", version=VERSION)
    parser.add_argument("-o", "--output-dir", dest="output_report_dir",
                        default=None)
    args = parser.parse_args()
    if args.output_report_dir is None:
        args.output_report_dir = os.path.join(os.getcwd(), "camsa_report_{date}".format(
            date=datetime.datetime.now().strftime("%b_%d_%Y__%H_%M")))

    aps = defaultdict(list)
    for file_name in args.input:
        file_name = os.path.abspath(os.path.expanduser(file_name))
        with open(file_name, "rt") as source:
            camsa_io.read_pairs(source=source, delimiter="\t", destination=aps,
                                default_cw_eae=args.cw_exact, default_cw_pae=args.cw_prob)

    or_ap_by_id = {}
    merged_ap_by_id = {}
    merged_aps = merge_aps(aps=aps, or_ap_by_id=or_ap_by_id, m_ap_by_id=merged_ap_by_id)

    aps_source = defaultdict(list)
    for ap in merged_aps:
        if ap.cw < args.min_cw:
            continue
        for u, v in ap.get_edges():
            aps_source[u].append(ap)
            aps_source[v].append(ap)

    tmp_individual_assemblies = defaultdict(list)
    for ap in merged_aps:
        for source_name in ap.sources:
            tmp_individual_assemblies[source_name].append(ap)

    individual_assemblies = [Assembly(name=name, aps=aps) for name, aps in tmp_individual_assemblies.items()]
    grouped_assemblies = []
    for i in range(1, len(individual_assemblies) + 1):
        for assembly_combination in itertools.combinations(sorted([a.name for a in individual_assemblies]), i):
            assembly = Assembly(name=assembly_combination, aps=[])
            for ap in merged_aps:
                if len(ap.sources) == len(assembly_combination) and all(map(lambda entry: entry in ap.sources, assembly_combination)):
                    assembly.aps.append(ap)
            grouped_assemblies.append(assembly)
    grouped_assemblies = list(filter(lambda a: len(a.aps) > 0, sorted(grouped_assemblies, key=lambda entry: len(entry.aps), reverse=True)))

    ag = AssemblyGraph()
    for ap in merged_aps:
        for (u, v) in ap.get_edges():
            ag.add_edge(u=u, v=v, weight=ap.cw)
    max_non_conflicting_assembly_graph = ag.get_maximal_non_conflicting_assembly_graph()

    for ap in merged_aps:
        participates = False
        for u, v in ap.get_edges():
            if max_non_conflicting_assembly_graph.has_edge(u=u, v=v):
                participates = True
                ap.participation_ctg1_or = "+" if u.endswith("h") else "-"
                ap.participation_ctg2_or = "+" if v.endswith("t") else "-"
                break
        ap.participates_in_max_non_conflicting_assembly = participates
        if not participates:
            conflicted_aps, semi_conflicted_aps = get_all_conflicted_aps(ap=ap, all_aps=aps_source)
            for c_ap in conflicted_aps:
                inter = set(ap.sources).intersection(set(c_ap.sources))
                for source in inter:
                    if source not in c_ap.in_conflicted:
                        c_ap.in_conflicted.append(source)
                    if source not in ap.in_conflicted:
                        ap.in_conflicted.append(source)
                other_out_c_sources = set(ap.sources).difference(c_ap.sources)
                for source in other_out_c_sources:
                    if source not in c_ap.out_conflicted:
                        c_ap.out_conflicted.append(source)
                self_out_c_sources = set(c_ap.sources).difference(ap.sources)
                for source in self_out_c_sources:
                    if source not in ap.out_conflicted:
                        ap.out_conflicted.append(source)
            for c_ap in semi_conflicted_aps:
                inter = set(ap.sources).intersection(set(c_ap.sources))
                for source in inter:
                    if source not in c_ap.in_semi_conflicted:
                        c_ap.in_semi_conflicted.append(source)
                    if source not in ap.in_semi_conflicted:
                        ap.in_semi_conflicted.append(source)
                other_out_c_sources = set(ap.sources).difference(c_ap.sources)
                for source in other_out_c_sources:
                    if source not in c_ap.out_semi_conflicted:
                        c_ap.out_semi_conflicted.append(source)
                self_c_sources = set(ap.sources).difference(c_ap.sources)
                for source in self_c_sources:
                    if source not in ap.out_semi_conflicted:
                        ap.out_semi_conflicted.append(source)

    args.output_report_dir = os.path.abspath(os.path.expanduser(args.output_report_dir))
    if not os.path.exists(args.output_report_dir):
        os.makedirs(args.output_report_dir)

    if os.path.exists(os.path.join(args.output_report_dir, "libs")) and os.path.isdir(
            os.path.join(args.output_report_dir, "libs")):
        shutil.rmtree(os.path.join(args.output_report_dir, "libs"))
    root_dir = os.path.dirname(os.path.realpath(__file__))
    shutil.copytree(src=os.path.join(root_dir, "libs"), dst=os.path.join(args.output_report_dir, "libs"))
    output_html_report_file_name = os.path.join(args.output_report_dir, args.output_report_file)

    template = Template(source=open(os.path.join(root_dir, "report_template.html"), "rt").read())

    individual_assemblies.sort(key=lambda assembly: assembly.name.lower())
    source_names_to_ids = {assembly.name: cnt for cnt, assembly in enumerate(individual_assemblies, start=1)}
    sources = [assembly.name for assembly in individual_assemblies]

    source_colors = {source_names_to_ids[source]: color for source, color in zip(sources, ['red', 'blue', 'green', 'purple',
                                                                                           'orange', 'pink', 'brown', 'navy', 'steelblue'])}

    with open(output_html_report_file_name, "wt") as dest:
        print(template.render(
            data={
                "assemblies": individual_assemblies,
                "assemblies_intersections": [],
                "assemblies_conflicts": [],
                "graph_compiled": False,
                "aps": merged_aps,
                "names_to_id": source_names_to_ids,
                "source_colors": source_colors,
                "grouped_assemblies": grouped_assemblies
            },
            meta={
                "camsa": {
                    "version": VERSION
                }
            }),
            file=dest)
