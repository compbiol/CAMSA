#! /usr/bin/env python
# -*- coding: utf-8 -*-
import datetime
import logging
import os
import sys
from collections import defaultdict

import configargparse

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))))

import camsa
import camsa.utils.ragout.io as ragout_io
from camsa.core.data_structures import AssemblyPoint
import camsa.core.io as camsa_io
from camsa.utils.ragout.shared import filter_indels, filter_duplications
from camsa.utils.ragout.shared import filter_blocks_by_good_genomes, filter_blocks_by_bad_genomes, get_all_genomes_from_blocks


def get_assembly_points(seq_of_blocks):
    result = []
    for l_block, r_block in zip(seq_of_blocks[:-1], seq_of_blocks[1:]):
        ap = AssemblyPoint(seq1=l_block.name, seq2=r_block.name, seq1_or=l_block.strand, seq2_or=r_block.strand, sources={l_block.parent_seq.genome_name, r_block.parent_seq.genome_name}, gap_size=r_block.start - l_block.end)
        result.append(ap)
    return result


if __name__ == "__main__":
    full_description = camsa.full_description_template.format(
        names=camsa.CAMSA_AUTHORS,
        affiliations=camsa.AFFILIATIONS,
        dummy=" ",
        tool="Converting Ragout formatted blocks into CAMSA assembly points.",
        information="For more information refer to {docs}".format(docs=camsa.CAMSA_DOCS_URL),
        contact=camsa.CONTACT)
    full_description = "=" * 80 + "\n" + full_description + "=" * 80 + "\n"
    parser = configargparse.ArgParser(description=full_description, formatter_class=configargparse.RawTextHelpFormatter,
                                      default_config_files=[os.path.join(camsa.root_dir, "logging.ini"),
                                                            os.path.join(camsa.root_dir, "utils", "ragout", "ragout_coords2camsa_points.ini")])
    parser.add_argument("-c", "--config", is_config_file=True, help="Config file overwriting some of the default settings as well as any flag starting with \"--\".")

    parser.add_argument("--version", action="version", version=camsa.VERSION)
    parser.add_argument("ragout_coords", type=str, help="A path to ragout coords file")
    parser.add_argument("--filter-indels", action="store_true", dest="filter_indels", default=False,
                        help="A flag to indicate whether to filter out all synteny blocks, that are not present exactly once in each of the good (all-bad) genomes.")
    parser.add_argument("--no-filter-indels", action="store_false", dest="filter_indels", default=False,
                        help="A flag to indicate whether to filter out all synteny blocks, that are not present exactly once in each of the good (all-bad) genomes.")
    parser.add_argument("--filter-duplications", action="store_true", dest="filter_duplications", default=True,
                        help="A flag to indicate whether to filter out all synteny blocks, that are present more than once in at least on of the good (all-bad) genomes.")
    parser.add_argument("--no-filter-duplications", action="store_false", dest="filter_duplications", default=True,
                        help="A flag to indicate whether to filter out all synteny blocks, that are present more than once in at least on of the good (all-bad) genomes.")
    parser.add_argument("--good-genomes", type=str, default="", help="A coma separated list of genome names, to be processed and conversed.\nDEFAULT: \"\" (i.e., all genomes are good)")
    parser.add_argument("--bad-genomes", type=str, default="Anc0", help="A coma separated list of genome names, to be excluded from processing and conversion.\nDEFAULT: \"\" (i.e., no genomes are bad)")
    parser.add_argument("-o", "--output", type=configargparse.FileType("wt"), default=sys.stdout)
    parser.add_argument("--o-genomes", type=str, dest="output_genomes", default="",
                        help="A coma separated list of genome names, which will determine the order inferred assembly points to be output.\nDEFAULT: \"\" (i.e., sorted list of good (all - bad) genomes)")
    parser.add_argument("--o-format", type=str, help="")
    parser.add_argument("--o-delimiter", default="\t", type=str, help="")
    parser.add_argument("--c-logging-level", dest="c_logging_level", default=logging.INFO, type=int,
                        choices=[logging.NOTSET, logging.DEBUG, logging.INFO, logging.WARNING, logging.ERROR, logging.CRITICAL],
                        help="Logging level for the converter.\nDEFAULT: {info}".format(info=logging.INFO))
    parser.add_argument("--c-logging-formatter-entry",
                        help="Format string for python logger.")
    args = parser.parse_args()

    start_time = datetime.datetime.now()

    logger = logging.getLogger("CAMSA.utils.ragout_coords2camsa_points")
    ch = logging.StreamHandler()
    ch.setLevel(args.c_logging_level)
    logger.setLevel(args.c_logging_level)
    logger.addHandler(ch)
    logger.info(full_description)
    logger.info(parser.format_values())
    ch.setFormatter(logging.Formatter(args.c_logging_formatter_entry))
    logger.info("Starting the converting process")

    sequences_by_ids, blocks_by_ids = ragout_io.read_from_file(path=args.ragout_coords, silent_fail=False, delimiter="\t")
    all_genomes = get_all_genomes_from_blocks(blocks_as_ids=blocks_by_ids)
    if args.good_genomes != "":
        args.good_genomes = set(args.good_genomes.split(","))
        filter_blocks_by_good_genomes(blocks_by_ids=blocks_by_ids, good_genomes=args.good_genomes)
    if args.bad_genomes != "":
        args.bad_genomes = set(args.bad_genomes.split(","))
        filter_blocks_by_bad_genomes(blocks_by_ids=blocks_by_ids, bad_genomes=args.bad_genomes)
    if args.filter_indels:
        filter_indels(blocks_by_ids=blocks_by_ids, all_genomes_as_set=(all_genomes if len(args.good_genomes) == 0 else args.good_genomes) - args.bad_genomes)
    if args.filter_duplications:
        filter_duplications(blocks_by_ids=blocks_by_ids)
    all_filtered_genomes = get_all_genomes_from_blocks(blocks_as_ids=blocks_by_ids)
    if args.output_genomes == "":
        output_genomes = sorted(all_filtered_genomes)
    else:
        output_genomes = args.output_genomes.split(",")
    for genome in output_genomes:
        if genome not in all_filtered_genomes:
            logger.critical("Genome {genome_name} specified with the --o-genomes flag was not present in all filtered genome {filtered_genomes}"
                            "".format(genome_name=genome, filtered_genomes=",".join(all_filtered_genomes)))
            exit(1)
    blocks_to_convert = []
    aps = []
    for genome in output_genomes:
        for block_id in sorted(blocks_by_ids.keys()):
            blocks = blocks_by_ids[block_id]
            blocks_by_genome = [block for block in blocks if block.parent_seq.genome_name == genome]
            blocks_to_convert.extend(blocks_by_genome)

        blocks_by_seq_ids = defaultdict(list)
        for block in blocks_to_convert:
            blocks_by_seq_ids[block.parent_seq.seq_name].append(block)

        for seq_id in list(blocks_by_seq_ids.keys()):
            blocks_by_seq_ids[seq_id] = sorted(blocks_by_seq_ids[seq_id], key=lambda block: (block.start, block.end))

        for seq_id in list(blocks_by_seq_ids.keys()):
            seq_of_blocks = blocks_by_seq_ids[seq_id]
            seq_aps = get_assembly_points(seq_of_blocks=seq_of_blocks)
            aps.extend(seq_aps)
    logger.info("Writing output to file \"{file_name}\"".format(file_name=args.output))
    camsa_io.write_assembly_points(assembly_points=aps, destination=args.output, output_setup=args.o_format, delimiter=args.o_delimiter)
    logger.info("Elapsed time: {el_time}".format(el_time=str(datetime.datetime.now() - start_time)))
