#! /usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function, division

import datetime
import logging
import numbers
import os
import sys
from collections import defaultdict, Counter

import Bio
import configargparse
import itertools
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
from camsa.utils.fasta.data_structures import IntraGapFilling, FlankingGapFilling


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


def get_fasta_for_meta_seq(meta_seq, orientation, fasta_by_ids):
    seq = fasta_by_ids[meta_seq.parent_seq_id].seq[meta_seq.start:meta_seq.end]
    seq = seq if meta_seq.strand == "+" else seq.reverse_complement()
    seq = seq if orientation == "+" else seq.reverse_complement()
    return seq


def get_meta_flanking_by_meta_seq_name(meta_seq_name, orientation, meta_sequences_by_name, meta_sequences_by_parent_ids, fasta_by_ids, extension_type):
    """

    :param meta_seq_name:
    :param orientation:
    :param meta_sequences_by_name:
    :param meta_sequences_by_parent_ids:
    :param fasta_by_ids:
    :param extension_type: "before" / "after"
    :return:
    """
    meta_sequences = meta_sequences_by_name[meta_seq_name]
    result = []
    for meta_sequence in meta_sequences:
        flanking_meta_seq = get_meta_flanking_by_meta_seq(meta_sequence=meta_sequence, extension_type=extension_type, fasta_by_ids=fasta_by_ids, meta_sequences_by_parent_ids=meta_sequences_by_parent_ids, orientation=orientation)
        if flanking_meta_seq is not None and flanking_meta_seq.length > 0:
            result.append(flanking_meta_seq)
    return result


def get_meta_flanking_by_meta_seq(extension_type, fasta_by_ids, meta_sequence, meta_sequences_by_parent_ids, orientation):
    if meta_sequence.name == meta_sequences_by_parent_ids[meta_sequence.parent_seq_id][0].name and meta_sequence.start > 0:
        if extension_type == "before" and meta_sequence.strand == orientation:
            extension_meta_seq = Sequence(name="extension",
                                          parent_seq_id=meta_sequence.parent_seq_id,
                                          start=0,
                                          end=meta_sequence.start,
                                          strand="+",
                                          seq_group_id=meta_sequence.seq_group_id)

            return extension_meta_seq
        elif extension_type == "after" and meta_sequence.strand != orientation:
            extension_meta_seq = Sequence(name="extension",
                                          parent_seq_id=meta_sequence.parent_seq_id,
                                          start=0,
                                          end=meta_sequence.start,
                                          strand="-",
                                          seq_group_id=meta_sequence.seq_group_id)

            return extension_meta_seq
    elif meta_sequence.name == meta_seqs_by_parent_ids[meta_sequence.parent_seq_id][-1].name and meta_sequence.end < len(fasta_by_ids[meta_sequence.parent_seq_id]):
        if extension_type == "before" and meta_sequence.strand != orientation:
            extension_meta_seq = Sequence(name="extension",
                                          parent_seq_id=meta_sequence.parent_seq_id,
                                          start=meta_sequence.end,
                                          end=len(fasta_by_ids[meta_sequence.parent_seq_id]),
                                          strand="-",
                                          seq_group_id=meta_sequence.seq_group_id)
            return extension_meta_seq
        elif extension_type == "after" and meta_sequence.strand == orientation:
            extension_meta_seq = Sequence(name="extension",
                                          parent_seq_id=meta_sequence.parent_seq_id,
                                          start=meta_sequence.end,
                                          end=len(fasta_by_ids[meta_sequence.parent_seq_id]),
                                          strand="+",
                                          seq_group_id=meta_sequence.seq_group_id)
            return extension_meta_seq


def collinearity(meta_seq1, meta_seq2, orientation1, orientation2):
    if meta_seq1.strand == orientation1 and meta_seq2.strand == orientation2 and meta_seq1.end <= meta_seq2.start:
        return 1
    elif meta_seq1.strand == reverse_or(orientaiton=orientation1) and meta_seq2.strand == reverse_or(orientaiton=orientation2) and meta_seq2.end <= meta_seq1.start:
        return -1
    return 0


def get_dummy_gap_filling(gap_size, sep, default_gap_size):
    if gap_size == "?":
        gap_size = default_gap_size
    result = Seq(sep * int(gap_size))
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
                        help="A stream of fasta formatted sequences of fragments, that participate in the scaffold assembly represented in form of CAMSA points")
    parser.add_argument("--points", type=configargparse.FileType("rt"), required=True,
                        help="A stream of CAMSA formatted assembly points, representing a scaffold assembly, that is converted into FASTA formatted sequences")
    parser.add_argument("--seqi", default="", type=str,
                        help="A stream of CAMSA formatted information about fragments, that participate in the reported assembly points.\nDEFAULT: dummy seqi records are created for each fasta record, if non provided.")
    parser.add_argument("--seqi-delimiter", default="\t", type=str,
                        help="A delimiter character for the file containing sequences' information.\nDEFAULT: \\t")
    parser.add_argument("--genomes-order", default="", type=str,
                        help="Order in which genomes sequences (i.e., seq_group_id column in *.seqi file) shall be for resulting sequence reconstruction.\nDEFAULT: sorted order of all seq_group_id column values. Seqi entries without such attribute are assigned the \"Default\" seq_group and used last.")
    parser.add_argument("--extend-ends", action="store_true", dest="extend_ends",
                        help="Flag specifying whether to extend the extremities in the reconstructed assembly, when possible, using flanking extremities in the sequences, from which fragments come.\nDEFAULT: True", default=True)
    parser.add_argument("--no-extend-ends", action="store_false", dest="extend_ends",
                        help="Flag specifying whether to extend the extremities in the reconstructed assembly, when possible, using flanking extremities in the sequences, from which fragments come.\nDEFAULT: True", default=True)
    parser.add_argument("--fill-gaps", action="store_true", dest="fill_gaps",
                        help="Flag specifying whether to fill gaps in the reconstructed assembly using information in the sequences, from which fragments participating in assembly points come.\nDEFAULT: True", default=True)
    parser.add_argument("--no-fill-gaps", action="store_false", dest="fill_gaps",
                        help="Flag specifying whether to fill gaps in the reconstructed assembly using information in the sequences, from which fragments participating in assembly points come.\nDEFAULT: True", default=True)
    parser.add_argument("--fill-gaps-unknown", action="store_true",
                        help="A flag specifying whether to fill gaps with unknown gap size or nor.\nDEFAULT: False (i.e., Ns will be inserted)")
    parser.add_argument("--gap-diff-threshold-per", default=10.0, type=float,
                        help="A percentage wise difference (i.e., in [0, 100]) specifying how much the length of the filling can differ from the given gap size.\nDEFAULT: 10")
    parser.add_argument("--gap-diff-threshold-bp", default=1000, type=float,
                        help="A bp wise difference (i.e., [0,..]) specifying how much the length of the filling can differ from the given gap size.\nDEFAULT: 10")
    parser.add_argument("--gap-diff-threshold-max", action="store_true", dest="gap_diff_max_over_min", default=False,
                        help="A flag specifying whether to take maximum possible difference (i.e., between percentage wise and bp wise) or a minimum one.\nDEFAULT: False (i.e., minimum, more conservative, difference is taken as a threshold)")
    parser.add_argument("--gap-fill-no-intra-first", action="store_false", dest="gap_fill_intra_first", default=True,
                        help="A flag specifying whether gaps shall be fist filled with intra-fragment sequence, rather than using two flanking, possible overlapping sequences.\nDEFAULT: True")
    parser.add_argument("--gap-fill-no-inter", action="store_false", dest="gap_fill_inter", default=True,
                        help="A flag specifying whether it is allowed to fill gaps with pairs flanking sequences or not.")
    # parser.add_argument("--allow-singletons", action="store_true", dest="allow_singletons", default=False,
    #                     help="Whether to include scaffolds, that were not mentioned in the CAMSA formatted assembly points or in there \nDEFAULT: False")
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

    logger.info("Reading fasta of contigs/scaffolds from {file}".format(file=args.fasta))
    frag_fasta_by_id = {}
    s_cnt = 0
    for record in SeqIO.parse(args.fasta, "fasta"):
        frag_fasta_by_id[record.id] = record
        s_cnt += 1
    logger.debug("Obtained {s_cnt} fasta fragments".format(s_cnt=s_cnt))

    sequences_by_ids = defaultdict(list)
    if args.seqi != "":
        logger.debug("Obtaining \"meta\"-sequence info for the fragments, that can be involved in the assembly points, from {file}".format(file=args.seqi))
        with open(args.seqi, "rt") as source:
            read_seqi_from_input_sources(source=source, delimiter=args.seqi_delimiter, destination=sequences_by_ids)
        logger.debug("Obtained {ms_cnt} meta-sequences records".format(ms_cnt=len(sequences_by_ids)))

    participating_sequences_by_ids = {}
    # for every fragment name that is not present as meta-sequence info a dumb one is created
    for fragment_aps in fragments:
        for f1, f1_or, f2, f2_or, gap_size in fragment_aps:
            if f1 not in sequences_by_ids:
                logger.debug("Fragment {f1} didn't have a meta-sequence information for it. Creating a dummy one.".format(f1=f1))
                sequences_by_ids[f1].append(Sequence(name=f1))
            if f2 not in sequences_by_ids:
                sequences_by_ids[f2].append(Sequence(name=f2))
                logger.debug("Fragment {f2} didn't have a meta-sequence information for it. Creating a dummy one.".format(f2=f2))
            participating_sequences_by_ids[f1] = sequences_by_ids[f1]
            participating_sequences_by_ids[f2] = sequences_by_ids[f2]

    all_seq_groups = set()
    for seq_id in list(participating_sequences_by_ids.keys()):
        for seq in participating_sequences_by_ids[seq_id]:
            all_seq_groups.add(seq.seq_group_id)
    if args.genomes_order == "":
        args.genomes_order = sorted(all_seq_groups, reverse=True)
    else:
        args.genomes_order = args.genomes_order.split(",")
    for genome in args.genomes_order:
        if genome not in all_seq_groups:
            logger.critical("Genome with name {genome} specified in the --genomes-order flag value was not present in the inferred set of sequences {seq_groups}"
                            "".format(genome=genome, seq_groups=",".join(sorted(all_seq_groups))))
            exit(1)
    # sorting participating sequences into "batches" denoted by seq_group_id field, with a default "None" value to it
    #   initial order of sequences within each batch is preserved w.r.t. the way it is listed in the *.seqi file
    for seq_id in list(participating_sequences_by_ids.keys()):
        in_genomes_order = [seq for seq in participating_sequences_by_ids[seq_id] if seq.seq_group_id in args.genomes_order]
        attributed_sequences = [seq for seq in participating_sequences_by_ids[seq_id] if seq.seq_group_id != "Default" and seq.seq_group_id not in args.genomes_order]
        non_attributed_sequences = [seq for seq in participating_sequences_by_ids[seq_id] if seq.seq_group_id == "Default"]
        participating_sequences_by_ids[seq_id] = []
        for genome in args.genomes_order:
            participating_sequences_by_ids[seq_id].extend([seq for seq in in_genomes_order if seq.seq_group_id == genome])
        attributed_sequences = sorted(attributed_sequences, reverse=True, key=lambda entry: entry.seq_group_id)
        participating_sequences_by_ids[seq_id].extend(attributed_sequences)
        participating_sequences_by_ids[seq_id].extend(non_attributed_sequences)

    meta_seqs_by_parent_ids = defaultdict(list)
    for seq_batch in participating_sequences_by_ids.values():
        for seq in seq_batch:
            meta_seqs_by_parent_ids[seq.parent_seq_id].append(seq)

    # for each parent sequence id all its "child" members are sorted with respect to their location on the parent sequence
    for parent_seq_id in list(meta_seqs_by_parent_ids.keys()):
        meta_seqs_by_parent_ids[parent_seq_id] = sorted(meta_seqs_by_parent_ids[parent_seq_id], key=lambda seq: (seq.start, seq.end))

    logger.info("Total number of fasta sequences is {seq_cnt}".format(seq_cnt=len(frag_fasta_by_id)))

    for seq_batch in participating_sequences_by_ids.values():
        for seq in seq_batch:
            if seq.parent_seq_id not in frag_fasta_by_id:
                logger.critical("Fragment {sequence} (parent for {block_name}) is not present in supplied fasta file. Exiting.".format(sequence=seq.parent_seq_id, block_name=seq.name))
                exit(1)

    used_fragments = set()
    filled_gaps_by_origin_cnt = Counter()
    extended_extremities_by_origin_cnt = Counter()
    total_gaps_cnt = 0
    total_extremities_cnt = 0
    for s_cnt, fragment_aps in enumerate(fragments):
        current = Seq("")
        total_extremities_cnt += 2
        for f_cnt, (f1, f1_or, f2, f2_or, gap_size) in enumerate(fragment_aps):
            main_meta_seq1 = participating_sequences_by_ids[f1][0]
            main_meta_seq2 = participating_sequences_by_ids[f2][0]

            seq1 = get_fasta_for_meta_seq(meta_seq=main_meta_seq1, orientation=f1_or, fasta_by_ids=frag_fasta_by_id)

            if f_cnt == 0 and args.extend_ends:
                logger.debug("Trying to extend the end for the new scaffold, based on fragment {f1}".format(f1=main_meta_seq1.name))
                flanking_meta_seqs = get_meta_flanking_by_meta_seq_name(meta_seq_name=main_meta_seq1.name, orientation=f1_or,
                                                                        meta_sequences_by_name=participating_sequences_by_ids, meta_sequences_by_parent_ids=meta_seqs_by_parent_ids,
                                                                        fasta_by_ids=frag_fasta_by_id, extension_type="before")
                if len(flanking_meta_seqs) > 0:
                    extension_meta_seq = flanking_meta_seqs[0]
                    extension = get_fasta_for_meta_seq(meta_seq=extension_meta_seq, orientation="+", fasta_by_ids=frag_fasta_by_id)
                    current += extension
                    extended_extremities_by_origin_cnt[extension_meta_seq.seq_group_id] += 1
                    logger.debug("Extended the end of the new scaffold based on the fragment {f1}, which appears to be the first in the \"chain\"".format(f1=f1))
                    logger.debug("\textension of length {ex_length} was obtained from the parent sequence {p_seq} in sequence group {seq_group_id}".format(ex_length=extension_meta_seq.length, p_seq=extension_meta_seq.parent_seq_id, seq_group_id=extension_meta_seq.seq_group_id))
                else:
                    logger.debug("Were unable to extend new scaffold based on the fragment {f1}, which appears to be the first in the \"chain\"".format(f1=f1))

            current += seq1
            total_gaps_cnt += 1
            filled_gap = False
            if args.fill_gaps:
                ap_gap_size = gap_size if isinstance(gap_size, numbers.Number) else "?"
                if ap_gap_size == "?" and not args.fill_gaps_unknown:
                    gap_filling = get_dummy_gap_filling(gap_size=ap_gap_size, sep=args.c_sep, default_gap_size=args.c_sep_length)
                    current += gap_filling
                    continue
                logger.debug("Trying to fill the gap for the new scaffold between fragments {f1} and {f2}".format(f1=f1, f2=f2))
                seq1_by_batches = defaultdict(list)
                seq1_batches_names = []
                for seq in participating_sequences_by_ids[f1]:
                    seq1_by_batches[seq.seq_group_id].append(seq)
                    if seq.seq_group_id not in seq1_batches_names:
                        seq1_batches_names.append(seq.seq_group_id)
                seq2_by_batches = defaultdict(list)
                seq2_batches_names = []
                for seq in participating_sequences_by_ids[f2]:
                    seq2_by_batches[seq.seq_group_id].append(seq)
                    if seq.seq_group_id not in seq2_batches_names:
                        seq2_batches_names.append(seq.seq_group_id)
                common_batches = [seq_group_id for seq_group_id in seq1_batches_names if seq_group_id in seq2_batches_names]
                logger.debug("{batches_cnt} sequence groups were identified where both {f1} and {f2} have sequence information".format(batches_cnt=len(common_batches), f1=f1, f2=f2))

                best_inter_gap_filling = None
                for seq_group_id in common_batches:
                    logger.debug("Trying to fill the gap using sequence group {seq_group_id}".format(seq_group_id=seq_group_id))
                    meta_seq1s = seq1_by_batches[seq_group_id]
                    meta_seq2s = seq2_by_batches[seq_group_id]
                    possible_fillings = []
                    for meta_seq1, meta_seq2 in itertools.product(meta_seq1s, meta_seq2s):
                        # case when both sequences reported in the assembly point have flanking sequences on their parent fragments
                        meta_seq1_flanking = get_meta_flanking_by_meta_seq(extension_type="after", fasta_by_ids=frag_fasta_by_id, meta_sequence=meta_seq1, meta_sequences_by_parent_ids=meta_seqs_by_parent_ids, orientation=f1_or)
                        meta_seq2_flanking = get_meta_flanking_by_meta_seq(extension_type="before", fasta_by_ids=frag_fasta_by_id, meta_sequence=meta_seq2, meta_sequences_by_parent_ids=meta_seqs_by_parent_ids, orientation=f2_or)
                        if meta_seq1_flanking is None or meta_seq2_flanking is None:
                            pass
                        else:
                            start_seq = get_fasta_for_meta_seq(meta_seq=meta_seq1_flanking, orientation="+", fasta_by_ids=frag_fasta_by_id)
                            end_seq = get_fasta_for_meta_seq(meta_seq=meta_seq2_flanking, orientation="+", fasta_by_ids=frag_fasta_by_id)
                            gap_filling = FlankingGapFilling(seq1=start_seq, seq2=end_seq, seq=None)
                            if args.gap_fill_inter:
                                possible_fillings.append(gap_filling)

                        # case when two sequences reported in the assembly point are on the same parent assembly and are co-oriented w.r.t. relative orientation reported in the assembly point
                        if meta_seq1.parent_seq_id == meta_seq2.parent_seq_id:
                            formation = collinearity(meta_seq1=meta_seq1, meta_seq2=meta_seq2, orientation1=f1_or, orientation2=f2_or)
                            if formation == 0:
                                continue
                            if formation == 1:
                                meta_seq_filling = Sequence(name="filling", parent_seq_id=meta_seq1.parent_seq_id, start=meta_seq1.end, end=meta_seq2.start, strand="+", seq_group_id=seq_group_id)
                            else:
                                meta_seq_filling = Sequence(name="filling", parent_seq_id=meta_seq1.parent_seq_id, start=meta_seq2.end, end=meta_seq1.start, strand="-", seq_group_id=seq_group_id)
                            filling_seq = get_fasta_for_meta_seq(meta_seq=meta_seq_filling, orientation="+", fasta_by_ids=frag_fasta_by_id)
                            gap_filling = IntraGapFilling(seq=filling_seq)
                            possible_fillings.append(gap_filling)
                    if len(possible_fillings) == 0:
                        logger.debug("Failed to fill the gap using sequence group {seq_group_id}".format(seq_group_id=seq_group_id))
                    else:
                        suitable_gap_fillings = []
                        per_based_gap_error = (ap_gap_size / 100 * args.gap_diff_threshold_per) if ap_gap_size != "?" else 0
                        bp_based_gap_error = args.gap_diff_threshold_bp if ap_gap_size != "?" else 0
                        gap_error = max(per_based_gap_error, bp_based_gap_error) if args.gap_diff_max_over_min else min(per_based_gap_error, bp_based_gap_error)
                        if ap_gap_size != "?":
                            gap_error = min(ap_gap_size, gap_error)
                            logger.debug("Inferred gap size error is {gse} for the gap of size {gs}".format(gse=gap_error, gs=ap_gap_size))
                        for gap_filling in possible_fillings:
                            if isinstance(gap_filling, IntraGapFilling) and gap_filling.suitable_to_fill_the_gap(gap_size=ap_gap_size, gap_size_error=gap_error):
                                gap_filling.compute_score(gap_size=ap_gap_size, gap_size_error=gap_error)
                                suitable_gap_fillings.append(gap_filling)
                            elif isinstance(gap_filling, FlankingGapFilling) and gap_filling.suitable_to_fill_the_gap(gap_size=ap_gap_size, gap_size_error=gap_error):
                                gap_filling.prepare_seq(gap_size=ap_gap_size, gap_size_error=gap_error)
                                gap_filling.compute_score(gap_size=ap_gap_size, gap_size_error=gap_error)
                                suitable_gap_fillings.append(gap_filling)
                        if len(suitable_gap_fillings) == 0:
                            logger.debug("Failed to fill the gap using sequence group {seq_group_id}".format(seq_group_id=seq_group_id))
                        else:
                            logger.debug("Filled the gap between fragments {f1} and {f2} using sequence group {seq_group_id}".format(f1=f1, f2=f2, seq_group_id=seq_group_id))
                            if args.gap_fill_intra_first:
                                intra_gap_fillings = [filling for filling in suitable_gap_fillings if isinstance(filling, IntraGapFilling)]
                                if len(intra_gap_fillings) > 0:
                                    gap_filling = sorted(suitable_gap_fillings, key=lambda entry: entry.score)[-1]
                                    current += gap_filling.seq
                                    filled_gap = True
                                    filled_gaps_by_origin_cnt[seq_group_id] += 1
                                else:
                                    if best_inter_gap_filling is None:
                                        best_inter_gap_filling = sorted(suitable_gap_fillings, key=lambda entry: entry.score)[-1]
                            else:
                                gap_filling = sorted(suitable_gap_fillings, key=lambda entry: entry.score)[-1]
                                current += gap_filling.seq
                                filled_gap = True
                                filled_gaps_by_origin_cnt[seq_group_id] += 1
                    if filled_gap:
                        break
                # is only triggered if the gap was NOT successfully filled with a non "N" sequence
                else:
                    if best_inter_gap_filling is not None:
                        gap_filling = best_inter_gap_filling.seq
                    else:
                        gap_filling = get_dummy_gap_filling(gap_size=ap_gap_size, sep=args.c_sep, default_gap_size=args.c_sep_length)
                    current += gap_filling

            if f_cnt == len(fragment_aps) - 1:
                seq2 = get_fasta_for_meta_seq(meta_seq=main_meta_seq2, orientation=f2_or, fasta_by_ids=frag_fasta_by_id)

                current += seq2

                if args.extend_ends:
                    logger.debug("Trying to extend the end for the new scaffold, based on fragment {f2}".format(f2=f2))
                    flanking_meta_seqs = get_meta_flanking_by_meta_seq_name(meta_seq_name=main_meta_seq2.name, orientation=f2_or,
                                                                            meta_sequences_by_name=participating_sequences_by_ids, meta_sequences_by_parent_ids=meta_seqs_by_parent_ids,
                                                                            fasta_by_ids=frag_fasta_by_id, extension_type="after")
                    if len(flanking_meta_seqs) > 0:
                        extension_meta_seq = flanking_meta_seqs[0]
                        extension = get_fasta_for_meta_seq(meta_seq=extension_meta_seq, orientation="+", fasta_by_ids=frag_fasta_by_id)
                        current += extension
                        extended_extremities_by_origin_cnt[extension_meta_seq.seq_group_id] += 1
                        logger.debug("Extended the end of the new scaffold based on the fragment {f2}, which appears to be the last in the \"chain\"".format(f2=f2))
                        logger.debug("\textension of length {ex_length} was obtained from the parent sequence {p_seq} in sequence group {seq_group_id}".format(ex_length=extension_meta_seq.length, p_seq=extension_meta_seq.parent_seq_id, seq_group_id=extension_meta_seq.seq_group_id))
                    else:
                        logger.debug("Were unable to extend new scaffold based on the fragment {f2}, which appears to be the first in the \"chain\"".format(f2=f2))

        name = args.scaffold_name_template.format(cnt=s_cnt)
        seq_record = SeqRecord(seq=current, id=name, description="")
        SeqIO.write(sequences=seq_record, handle=args.output, format="fasta")
    logger.info("Filled {gap_fill_cnt} / {total_gaps_cnt} gaps.".format(gap_fill_cnt=sum(filled_gaps_by_origin_cnt.values()), total_gaps_cnt=total_gaps_cnt))
    for key, value in filled_gaps_by_origin_cnt.items():
        logger.debug("\tFilled {gap_filled_cnt} gaps with source \"{seq_group_id}\"".format(gap_filled_cnt=value, seq_group_id=key))
    logger.info("Extended {end_ext_cnt} / {total_extremities_cnt} extremities.".format(end_ext_cnt=sum(extended_extremities_by_origin_cnt.values()), total_extremities_cnt=total_extremities_cnt))
    for key, value in extended_extremities_by_origin_cnt.items():
        logger.debug("\tExtended {extended_extremities_cnt} extremities of new scaffolds with source \"{seq_group_id}\"".format(extended_extremities_cnt=value, seq_group_id=key))
    # if args.allow_singletons:
    #     logger.info("Adding singleton fragments, that did not participate in any assembly points to the resulting assembly")
    #     for f_id, fragment in frag_fasta_by_id.items():
    #         if f_id not in used_fragments:
    #             SeqIO.write(sequences=fragment, handle=args.output, format="fasta")
    logger.info("New scaffolds were written to {file_name}".format(file_name=args.output))
    logger.info("All done!")
    logger.info("Elapsed time: {el_time}".format(el_time=str(datetime.datetime.now() - start_time)))
