# -*- coding: utf-8 -*-
import networkx as nx
from camsa.core.data_structures import ScaffoldAssemblyGraph


def analyze_and_update_assembly_points_based_on_reference(assembly, ref_assembly):
    reference_graph = ScaffoldAssemblyGraph.from_assembly_points(assembly_points=ref_assembly.aps, reference=True)
    get_correct_assembly_points_and_update_them(assembly=assembly, ref_assembly_graph=reference_graph)
    ref_all_pairs_paths = nx.all_pairs_shortest_path(G=reference_graph.graph)


def get_correct_assembly_points_and_update_them(assembly, ref_assembly_graph):
    for ap in assembly.aps:
        for u, v in ap.get_edges(sort=True):
            if ref_assembly_graph.graph.has_edge(u=u, v=v):
                ap.ref_metrics.present = True
                ref_ap_ids = set()
                for data in ref_assembly_graph.graph[u][v].values():
                    if "ap_id" in data:
                        ref_ap_ids.add(data["ap_id"])
                ap.present_ref_ids = sorted(ref_ap_ids)