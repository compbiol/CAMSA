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
from camsa.utils.ragout.shared import filter_indels, filter_duplications
from camsa.utils.ragout.shared import filter_blocks_by_good_genomes, filter_blocks_by_bad_genomes, get_all_genomes_from_blocks

if __name__ == "__main__":
    full_description = camsa.full_description_template.format(
        names=camsa.CAMSA_AUTHORS,
        affiliations=camsa.AFFILIATIONS,
        dummy=" ",
        tool="Computing coverage report for Ragout blocks.",
        information="For more information refer to {docs}".format(docs=camsa.CAMSA_DOCS_URL),
        contact=camsa.CONTACT)
    full_description = "=" * 80 + "\n" + full_description + "=" * 80 + "\n"
    parser = configargparse.ArgParser(description=full_description, formatter_class=configargparse.RawTextHelpFormatter,
                                      default_config_files=[os.path.join(camsa.root_dir, "logging.ini")])
    parser.add_argument("-c", "--config", is_config_file=True, help="Config file overwriting some of the default settings as well as any flag starting with \"--\".")

    parser.add_argument("--version", action="version", version=camsa.VERSION)
    parser.add_argument("ragout_coords", type=str, help="A path to ragout coords file")
    parser.add_argument("--filter-indels", action="store_true", dest="filter_indels", default=False)
    parser.add_argument("--no-fragment-stats", action="store_false", dest="fragment_stats", default=True)
    parser.add_argument("--no-genome-stats", action="store_false", dest="genome_stats", default=True)
    parser.add_argument("--filter-duplications", action="store_true", dest="filter_duplications", default=False)
    parser.add_argument("--good-genomes", type=str, default="", help="A coma separated list of genome names, to be processed and conversed.\nDEFAULT: \"\" (i.e., all genomes are good)")
    parser.add_argument("--bad-genomes", type=str, default="", help="A coma separated list of genome names, to be excluded from processing and conversion.\nDEFAULT: \"\" (i.e., no genomes are bad)")
    parser.add_argument("-o", "--output", type=configargparse.FileType("wt"), default=sys.stdout)
    parser.add_argument("--c-logging-level", dest="c_logging_level", default=logging.INFO, type=int,
                        choices=[logging.NOTSET, logging.DEBUG, logging.INFO, logging.WARNING, logging.ERROR, logging.CRITICAL],
                        help="Logging level for the converter.\nDEFAULT: {info}".format(info=logging.INFO))
    parser.add_argument("--c-logging-formatter-entry",
                        help="Format string for python logger.")
    args = parser.parse_args()

    start_time = datetime.datetime.now()

    logger = logging.getLogger("CAMSA.utils.ragout_coords2fasta")
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

    genomes = defaultdict(lambda: defaultdict(list))
    for block_list in blocks_by_ids.values():
        for block in block_list:
            genomes[block.parent_seq.genome_name][block.parent_seq.ragout_id].append(block)

    fragment_cov = {}
    if args.fragment_stats:
        for genome_name in genomes.keys():
            for seq_id in genomes[genome_name]:
                seq = sequences_by_ids[seq_id]
                cumulative_blocks_length = sum(block.length for block in genomes[genome_name][seq_id])
                fragment_cov[seq_id] = cumulative_blocks_length * 100.0 / seq.length

    genome_cov = {}
    if args.genome_stats:
        for genome_name in genomes.keys():
            total_genome_length = 0
            total_blocks_length = 0
            for seq_id in genomes[genome_name]:
                seq = sequences_by_ids[seq_id]
                total_genome_length += seq.length
                total_blocks_length += sum(block.length for block in genomes[genome_name][seq_id])
            genome_cov[genome_name] = total_blocks_length * 100.0 / total_genome_length

    if args.genome_stats:
        print("-" * 80, file=args.output)
        for genome_name in sorted(genomes.keys()):
            print("For genome \"{genome_name}\" {cov:.2f}% of its length is covered by filtered blocks".format(genome_name=genome_name, cov=genome_cov[genome_name]), file=args.output)
    if args.fragment_stats:
        for genome_name in sorted(genomes.keys()):
            print("-"*80, file=args.output)
            print("Detailed coverage stats for fragments in genome \"{genome_name}\"".format(genome_name=genome_name), file=args.output)
            for seq_id in sorted(genomes[genome_name].keys()):
                print("For fragment \"{fragment_name}\" {cov:.2f}% of its length is covered by filtered blocks".format(fragment_name=sequences_by_ids[seq_id].seq_name, cov=fragment_cov[seq_id]))

    logger.info("All done!")
    end_time = datetime.datetime.now()
    logger.info("Elapsed time: {el_time}".format(el_time=str(end_time - start_time)))