#! /usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function, division

import datetime
import logging
import os
import sys
from collections import defaultdict

import Bio
import configargparse
import networkx
from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))))

import camsa
from camsa.core.io import read_pairs
from camsa.core.merging import get_scaffold_edges


def get_scaffold_name_from_vertex(v):
    return v[:-1]


def get_assembly_edge(graph):
    for edge in graph.edges():
        v1, v2 = edge
        if v1[:-1] != v2[:-1]:
            return edge
    return None, None


def get_sequence_of_fragments_from_path(path):
    path, path_type = path
    if len(path) == 1:
        logger.error("A sequence, that contains a half of a scaffold. Something went wrong.")
        exit(1)
    result = []
    for scaf_extremity_vertex in path[:2]:
        orientation = "+" if scaf_extremity_vertex.endswith("t") else "-"
        result.append((get_scaffold_name_from_vertex(v=scaf_extremity_vertex), orientation))
    return result


if __name__ == "__main__":
    full_description = camsa.full_description_template.format(
        names=camsa.CAMSA_AUTHORS,
        affiliations=camsa.AFFILIATIONS,
        dummy=" ",
        tool="Converting CAMSA formatted scaffolding results into FASTA files.",
        information="For more information refer to wiki at github.com/aganezov/camsa/wiki",
        contact=camsa.CONTACT)
    full_description = "=" * 80 + "\n" + full_description + "=" * 80 + "\n"
    parser = configargparse.ArgParser(description=full_description, formatter_class=configargparse.RawTextHelpFormatter,
                                      default_config_files=[os.path.join(camsa.root_dir, "utils", "fasta", "camsa_points2fasta.ini"),
                                                            os.path.join(camsa.root_dir, "logging.ini")])
    parser.add_argument("-c", "--config", is_config_file=True, help="Config file overwriting some of the default settings as well as any flag starting with \"--\".")

    parser.add_argument("--fasta", type=configargparse.FileType("rt"), required=True)
    parser.add_argument("--points", type=configargparse.FileType("rt"), required=True)
    parser.add_argument("--c-sep", type=str)
    parser.add_argument("--c-sep-length", type=int)
    parser.add_argument("--scaffold-name-template", type=str)
    parser.add_argument("-o", "--output", type=configargparse.FileType("wt"), default=sys.stdout)

    parser.add_argument("--c-logging-level", dest="logging_level", default=logging.INFO, type=int,
                        choices=[logging.NOTSET, logging.DEBUG, logging.INFO, logging.WARNING, logging.ERROR, logging.CRITICAL],
                        help="Logging level for the converter.\nDEFAULT: {info}".format(info=logging.INFO))
    parser.add_argument("--c-logging-formatter-entry",
                        help="Format string for python logger.")
    args = parser.parse_args()
    start_time = datetime.datetime.now()

    logger = logging.getLogger("CAMSA.utils.camsa_points2fasta")
    ch = logging.StreamHandler()
    ch.setLevel(args.logging_level)
    logger.setLevel(args.logging_level)
    logger.addHandler(ch)
    logger.info(full_description)
    logger.info(parser.format_values())
    ch.setFormatter(logging.Formatter(args.c_logging_formatter_entry))
    logger.info("Starting the converting process")

    logger.info("Reading assembly points")
    assembly_points_by_sources = defaultdict(list)
    read_pairs(source=args.points, destination=assembly_points_by_sources)
    assembly_points_by_sources = [ap for ap_list in assembly_points_by_sources.values() for ap in ap_list]
    logger.info("A total of {ap_cnt} assembly points was obtained".format(ap_cnt=len(assembly_points_by_sources)))

    assembly_graph = networkx.Graph()

    scaffold_edges = get_scaffold_edges(assembly_points=assembly_points_by_sources)
    assembly_graph.add_edges_from(ebunch=scaffold_edges)

    for ap in assembly_points_by_sources:
        for (u, v, weight) in ap.get_edges(sort=True, weight=True):
            assembly_graph.add_edge(u, v)

    logger.debug("Checking that there are no in(semi)conflicting assembly points")
    for vertex in assembly_graph.nodes():
        degree = assembly_graph.degree(nbunch=vertex)
        if degree > 2:
            scaffold_name = get_scaffold_name_from_vertex(v=vertex)
            logger.error("Supplied assembly contained a conflict.")
            logger.error("Scaffold {scaffold_name} by its extremity {extremity_name} is reported as adjacent to more than one other scaffold's extremity"
                         "".format(scaffold_name=scaffold_name, extremity_name=vertex))
            exit(1)
    logger.debug("All clear, no (semi)conflicts, we can proceed")

    logger.info("Processing assemblies constructed from obtained assembly points")
    logger.debug("Extracting paths from assembly graph")
    paths = []
    for cc in networkx.connected_component_subgraphs(G=assembly_graph):
        origins = [v for v in cc.nodes() if cc.degree(v) == 1]
        if len(origins) == 2:
            path = networkx.shortest_path(G=cc, source=origins[0], target=origins[1])
            logger.debug("Extracted a linear scaffold of length {scaffold_length}, staring with {s_v} and ending with {e_v}"
                         "".format(scaffold_length=int(len(path) / 2),
                                   s_v=get_scaffold_name_from_vertex(v=origins[0]),
                                   e_v=get_scaffold_name_from_vertex(v=origins[1])))
            paths.append((path, "l"))
        if len(origins) == 1:
            logger.error("Something is wrong with the assembly graph. We have a connected component with a single vertex of degree 1.")
            exit(1)
        if len(origins) == 0:
            logger.debug("Encountered a circular chromosome. Splitting it at random assembly point")
            assembly_edge = get_assembly_edge(graph=cc)
            if assembly_edge[0] is None or assembly_edge[1] is None:
                logger.error("Something is wrong with the assembly graph. Couldn't find a scaffold edge in a circular scaffold.")
                exit(1)
            cc.remove_edge(u=assembly_edge[0], v=assembly_edge[1])
            path = networkx.shortest_path(G=cc, source=assembly_edge[0], target=assembly_edge[1])
            paths.append((path, "c"))
    logger.debug("Total number of extracted paths is {path_cnt}".format(path_cnt=len(paths)))
    logger.debug("Out of which {linear_cnt} are linear, and {circular_cnt} are circular"
                 "".format(linear_cnt=len([p for p in paths if p[1] == "l"]),
                           circular_cnt=len([p for p in paths if p[1] == "c"])))
    scaffolds = [get_sequence_of_fragments_from_path(path=p) for p in paths]
    logger.info("Total number of {scaffold_cnt} scaffolds was obtained from observed assembly points".format(scaffold_cnt=len(scaffolds)))

    logger.info("Reading fasta of contigs/scaffolds involved in the assembly points")
    frag_fasta_by_id = {}
    s_cnt = 0
    for record in SeqIO.parse(args.fasta, "fasta"):
        frag_fasta_by_id[record.id] = record
        s_cnt += 1
    logger.debug("Processed {cnt} records from fasta file \"{file_name}\"".format(cnt=s_cnt, file_name=args.fasta))
    logger.info("Total number of contig/scaffold sequences is {seq_cnt}".format(seq_cnt=len(frag_fasta_by_id)))

    for scaffold in scaffolds:
        for f, orientation in scaffold:
            if f not in frag_fasta_by_id:
                logging.critical("Fragment {f} which is present assembly points is not present in supplied fasta file. Exiting.")
                exit(1)

    logger.info("Outputting new scaffolds. Data is written to {file_name}".format(file_name=args.output))
    for s_cnt, scaffold in enumerate(scaffolds):
        current = Seq("")
        for f_cnt, (fragment, orientation) in enumerate(scaffold):
            if f_cnt > 0:
                current += Seq(args.c_sep * args.c_sep_length)
            if orientation == "+":
                current += frag_fasta_by_id[fragment].seq
            else:
                current += frag_fasta_by_id[fragment].reverse_complement().seq
        name = args.scaffold_name_template.format(cnt=s_cnt)
        seq_record = SeqRecord(seq=current, id=name, description="")
        SeqIO.write(sequences=seq_record, handle=args.output, format="fasta")
    logger.info("All done!")
    logger.info("Elapsed time: {el_time}".format(el_time=str(datetime.datetime.now() - start_time)))

