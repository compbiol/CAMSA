#! /usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function, division

import datetime
import logging
import numbers
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
from camsa.core.io import read_pairs, read_seqi_from_input_sources
from camsa.core.data_structures import get_scaffold_edges, Sequence


def get_scaffold_name_from_vertex(v):
    return v[:-1]


def reverse_or(orientaiton):
    return "+" if orientaiton == "-" else "-"


def reverse_ap(f1, f1_or, f2, f2_or):
    return f2, reverse_or(orientaiton=f2_or), f1, reverse_or(orientaiton=f1_or)


def get_assembly_edge(graph):
    for edge in graph.edges():
        v1, v2 = edge
        if v1[:-1] != v2[:-1]:
            return edge
    return None, None


def get_sequence_of_fragments_from_path(path, assembly_points_by_edges):
    path, path_type = path
    if len(path) < 2:
        logger.error("A sequence resembling an assembly points, that contains less than two scaffolds. Something went wrong.")
        exit(1)
    result = []
    for frag_extremity_v1, frag_extremity_v2 in zip(path[1::2], path[2::2]):
        f1_or = "+" if frag_extremity_v1.endswith("h") else "-"
        f2_or = "-" if frag_extremity_v2.endswith("h") else "+"
        ap = assembly_points_by_edges[tuple(sorted([frag_extremity_v1, frag_extremity_v2]))]
        gap_size = ap.gap_size
        result.append((get_scaffold_name_from_vertex(v=frag_extremity_v1), f1_or,
                       get_scaffold_name_from_vertex(v=frag_extremity_v2), f2_or,
                       gap_size))
    return result


if __name__ == "__main__":
    full_description = camsa.full_description_template.format(
        names=camsa.CAMSA_AUTHORS,
        affiliations=camsa.AFFILIATIONS,
        dummy=" ",
        tool="Converting CAMSA formatted scaffolding results into FASTA files.",
        information="For more information refer to {docs}".format(docs=camsa.CAMSA_DOCS_URL),
        contact=camsa.CONTACT)
    full_description = "=" * 80 + "\n" + full_description + "=" * 80 + "\n"
    parser = configargparse.ArgParser(description=full_description, formatter_class=configargparse.RawTextHelpFormatter,
                                      default_config_files=[os.path.join(camsa.root_dir, "utils", "fasta", "camsa_points2fasta.ini"),
                                                            os.path.join(camsa.root_dir, "logging.ini")])
    parser.add_argument("-c", "--config", is_config_file=True, help="Config file overwriting some of the default settings as well as any flag starting with \"--\".")

    parser.add_argument("--version", action="version", version=camsa.VERSION)
    parser.add_argument("--fasta", type=configargparse.FileType("rt"), required=True,
                        help="A stream of fasta formatted sequences of scaffolds, that participate in the scaffold assembly represented in form of CAMSA points")
    parser.add_argument("--points", type=configargparse.FileType("rt"), required=True,
                        help="A stream of CAMSA formatted assembly points, representing a scaffold assembly, that is converted into FASTA formatted sequences")
    parser.add_argument("--seqi", default="", type=str)
    parser.add_argument("--seqi-delimiter", default="\t", type=str)
    parser.add_argument("--fill-gaps", action="store_true")
    parser.add_argument("--extend-ends", action="store_true")
    parser.add_argument("--fill-gaps-unknown", action="store_true")
    parser.add_argument("--gap-by-ends-add", action="store_true")
    parser.add_argument("--gap-diff-threshold-per", default=10.0, type=float)
    parser.add_argument("--gap-diff-threshold-bp", default=1000, type=float)
    parser.add_argument("--allow-singletons", action="store_true", dest="allow_singletons", default=False,
                        help="Whether to include scaffolds, that were not mentioned in the CAMSA formatted assembly points\nDEFAULT: False")
    parser.add_argument("--c-sep", type=str,
                        help="A symbol, that is used to indicate gaps between scaffolds in the translated assemblies\nDEFAULT: N")
    parser.add_argument("--c-sep-length", type=int,
                        help="A default length, that is used for the gap size between scaffolds in the translated assemblies. Used in case, when gap-size column has \"?\" value\nDEFAULT: 20")
    parser.add_argument("--scaffold-name-template", type=str,
                        help="Python string template for the scaffold ids, in the produced FASTA formatted sequences. \"cnt\" attribute can be utilized\nDEFAULT: scaffold_{cnt}")
    parser.add_argument("-o", "--output", type=configargparse.FileType("wt"), default=sys.stdout,
                        help="A stream to which the FASTA formatted converted sequence, representing the CAMSA formatted scaffold assembly, is output\nDEFAULT: stdout")

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
    assembly_points_by_edges = {}

    for ap in assembly_points_by_sources:
        for (u, v, weight) in ap.get_edges(sort=True, weight=True):
            assembly_graph.add_edge(u, v)
            assembly_points_by_edges[tuple(sorted([u, v]))] = ap

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
    fragments = [get_sequence_of_fragments_from_path(path=p, assembly_points_by_edges=assembly_points_by_edges) for p in paths]
    logger.info("Total number of {scaffold_cnt} scaffolds was obtained from observed assembly points".format(scaffold_cnt=len(fragments)))

    logger.info("Reading fasta of contigs/scaffolds involved in the assembly points")
    frag_fasta_by_id = {}
    s_cnt = 0
    for record in SeqIO.parse(args.fasta, "fasta"):
        frag_fasta_by_id[record.id] = record
        s_cnt += 1
    sequences_by_ids = {}
    if args.seqi != "":
        with open(args.seqi, "rt") as source:
            read_seqi_from_input_sources(source=source, delimiter=args.seqi_delimiter, destination=sequences_by_ids)

    participating_sequences_by_ids = {}
    # for every fragment name that is not present as meta-sequence info a dumb one is created
    for fragment_aps in fragments:
        for f1, f1_or, f2, f2_or, gap_size in fragment_aps:
            if f1 not in sequences_by_ids:
                sequences_by_ids[f1] = Sequence(name=f1)
            if f2 not in sequences_by_ids:
                sequences_by_ids[f2] = Sequence(name=f2)
            participating_sequences_by_ids[f1] = sequences_by_ids[f1]
            participating_sequences_by_ids[f2] = sequences_by_ids[f2]

    meta_seqs_by_parent_ids = defaultdict(list)
    for seq in participating_sequences_by_ids.values():
        meta_seqs_by_parent_ids[seq.parent_seq_id].append(seq)

    for parent_seq_id in list(meta_seqs_by_parent_ids.keys()):
        meta_seqs_by_parent_ids[parent_seq_id] = sorted(meta_seqs_by_parent_ids[parent_seq_id], key=lambda seq: (seq.start, seq.end))

    logger.debug("Processed {cnt} records from fasta file \"{file_name}\"".format(cnt=s_cnt, file_name=args.fasta))
    logger.info("Total number of contig/scaffold sequences is {seq_cnt}".format(seq_cnt=len(frag_fasta_by_id)))

    for sequence in participating_sequences_by_ids.values():
        if sequence.parent_seq_id not in frag_fasta_by_id:
            logging.critical("Fragment {sequence} (parent for {block_name}) is not present in supplied fasta file. Exiting.".format(sequence=sequence.parent_seq_id, block_name=sequence.name))
            exit(1)

    used_fragments = set()
    logger.info("Outputting new scaffolds. Data is written to {file_name}".format(file_name=args.output))
    for s_cnt, fragment_aps in enumerate(fragments):
        current = Seq("")
        for f_cnt, (f1, f1_or, f2, f2_or, gap_size) in enumerate(fragment_aps):
            meta_seq1 = participating_sequences_by_ids[f1]
            meta_seq2 = participating_sequences_by_ids[f2]

            seq1 = frag_fasta_by_id[meta_seq1.parent_seq_id].seq[meta_seq1.start:meta_seq1.end]
            seq1 = seq1 if meta_seq1.strand == "+" else seq1.reverse_complement()
            seq1 = seq1 if f1_or == "+" else seq1.reverse_complement()

            if f_cnt == 0 and args.extend_ends:
                if meta_seq1.name == meta_seqs_by_parent_ids[meta_seq1.parent_seq_id][0].name \
                        and meta_seq1.strand == f1_or and meta_seq1.start > 0:
                    extension = frag_fasta_by_id[meta_seq1.parent_seq_id][:meta_seq1.start].seq
                    current += extension
                if meta_seq1.name == meta_seqs_by_parent_ids[meta_seq1.parent_seq_id][-1].name \
                        and meta_seq1.strand != f1_or:
                    extension = frag_fasta_by_id[meta_seq1.parent_seq_id][meta_seq1.end:].seq.reverse_complement()
                    current += extension

            current += seq1

            filled_gap = False
            if args.fill_gaps:
                if meta_seq1.parent_seq_id == meta_seq2.parent_seq_id:
                    if meta_seq1.end < meta_seq2.start and f1_or == meta_seq1.strand and f2_or == meta_seq2.strand:
                        gap = frag_fasta_by_id[meta_seq1.parent_seq_id][meta_seq1.end:meta_seq2.start].seq
                    elif meta_seq2.end < meta_seq1.start and f1_or == reverse_or(orientaiton=meta_seq1.strand) and f2_or == reverse_or(orientaiton=meta_seq2.strand):
                        gap = frag_fasta_by_id[meta_seq1.parent_seq_id][meta_seq2.end:meta_seq1.start].seq.reverse_complement()
                    else:
                        gap = Seq("")
                    gap_length = len(gap)
                    ap_gap_size = gap_size if isinstance(gap_size, numbers.Number) else "?"
                    if ap_gap_size == "?":
                        if args.fill_gaps_unknown:
                            current += gap
                            filled_gap = gap_length > 0
                    else:
                        diff = abs(ap_gap_size - gap_length)
                        if diff * 100.0 / ap_gap_size < args.gap_diff_threshold_per:
                            current += gap
                            filled_gap = gap_length > 0
                        elif diff < args.gap_diff_threshold_bp:
                            current += gap
                            filled_gap = gap_length > 0
                else:
                    f1_extension = Seq("")
                    f1_is_flanking = False
                    if meta_seq1.name == meta_seqs_by_parent_ids[meta_seq1.parent_seq_id][0].name \
                            and meta_seq1.strand != f1_or and meta_seq1.start > 0:
                        f1_extension = frag_fasta_by_id[meta_seq1.parent_seq_id][:meta_seq1.start].seq.reverse_complement()
                        f1_is_flanking = True
                    if meta_seq1.name == meta_seqs_by_parent_ids[meta_seq1.parent_seq_id][-1].name \
                            and meta_seq1.strand == f1_or:
                        f1_extension = frag_fasta_by_id[meta_seq1.parent_seq_id][meta_seq1.end:].seq
                        f1_is_flanking = True
                    f2_extension = Seq("")
                    f2_is_flanking = False
                    if meta_seq2.name == meta_seqs_by_parent_ids[meta_seq2.parent_seq_id][0].name \
                            and meta_seq2.strand == f2_or and meta_seq2.start > 0:
                        f2_extension = frag_fasta_by_id[meta_seq2.parent_seq_id].seq[:meta_seq2.start]
                        f2_is_flanking = True
                    if meta_seq2.name == meta_seqs_by_parent_ids[meta_seq2.parent_seq_id][-1].name \
                            and meta_seq2.strand != f2_or:
                        f2_extension = frag_fasta_by_id[meta_seq2.parent_seq_id].seq[meta_seq2.end:].reverse_complement()
                        f2_is_flanking = True

                    if f1_is_flanking and f2_is_flanking:
                        cumulative_length = len(f1_extension) + len(f2_extension)
                        ap_gap_size = gap_size if isinstance(gap_size, numbers.Number) else "?"
                        if ap_gap_size == "?":
                            if args.fill_gaps_unknown:
                                current += f1_extension + f2_extension
                                filled_gap = cumulative_length > 0
                        else:
                            if cumulative_length < ap_gap_size:
                                if args.gap_by_ends_add:
                                    remaining_length = ap_gap_size - cumulative_length
                                    current += f1_extension + Seq(args.c_sep*int(remaining_length)) + f2_extension
                                else:
                                    diff = abs(ap_gap_size - cumulative_length)
                                    if diff * 100.0 / ap_gap_size < args.gap_diff_threshold_per:
                                        current += f1_extension + f2_extension
                                        filled_gap = cumulative_length > 0
                                    elif diff < args.gap_diff_threshold_bp:
                                        current += f1_extension + f2_extension
                                        filled_gap = cumulative_length > 0
            if not filled_gap:
                gap_length = gap_size if isinstance(gap_size, numbers.Number) else args.c_sep_length
                if gap_length <= 0:
                    gap_length = args.c_sep_length
                gap = Seq(args.c_sep * int(gap_length))
                current += gap

            used_fragments.add(meta_seq1.parent_seq_id)
            used_fragments.add(meta_seq2.parent_seq_id)

            if f_cnt == len(fragment_aps) - 1:
                seq2 = frag_fasta_by_id[meta_seq2.parent_seq_id].seq[meta_seq2.start:meta_seq2.end]
                seq2 = seq2 if meta_seq2.strand == "+" else seq2.reverse_complement()
                seq2 = seq2 if f2_or == "+" else seq2.reverse_complement()
                current += seq2
                if args.extend_ends:
                    if meta_seq2.name == meta_seqs_by_parent_ids[meta_seq2.parent_seq_id][0].name \
                            and meta_seq2.strand != f2_or and meta_seq2.start > 0:
                        extension = frag_fasta_by_id[meta_seq2.parent_seq_id].seq[:meta_seq2.start].reverse_complement()
                        current += extension
                    if meta_seq2.name == meta_seqs_by_parent_ids[meta_seq2.parent_seq_id][-1].name \
                            and meta_seq2.strand == f2_or:
                        extension = frag_fasta_by_id[meta_seq2.parent_seq_id].seq[meta_seq2.end:]
                        current += extension

        name = args.scaffold_name_template.format(cnt=s_cnt)
        seq_record = SeqRecord(seq=current, id=name, description="")
        SeqIO.write(sequences=seq_record, handle=args.output, format="fasta")
    if args.allow_singletons:
        logger.info("Adding singleton fragments, that did not participate in any assembly points to the resulting assmebly")
        for f_id, fragment in frag_fasta_by_id.items():
            if f_id not in used_fragments:
                SeqIO.write(sequences=fragment, handle=args.output, format="fasta")
    logger.info("All done!")
    logger.info("Elapsed time: {el_time}".format(el_time=str(datetime.datetime.now() - start_time)))

