# -*- coding: utf-8 -*-
import networkx as nx
from camsa.core.data_structures import ScaffoldAssemblyGraph


def analyze_and_update_assembly_points_based_on_reference(assembly, ref_assembly):
    reference_graph = ScaffoldAssemblyGraph.from_assembly_points(assembly_points=ref_assembly.aps, reference=True)
    ref_all_pairs_paths = nx.all_pairs_shortest_path(G=reference_graph.graph)


def has_unoriented_assembly_points(assembly):
    return not all(map(lambda ap: ap.is_oriented, assembly.aps))



