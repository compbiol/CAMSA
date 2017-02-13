#! /usr/bin/env python
# -*- coding: utf-8 -*-
import datetime
import logging
import os
import sys

import configargparse

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))))

import camsa
from camsa.core.data_structures import Sequence
from camsa.core.io import write_seqi
from camsa.utils.ragout.shared import get_all_genomes_from_blocks, filter_blocks_by_good_genomes, filter_blocks_by_bad_genomes, filter_indels, filter_duplications
import camsa.utils.ragout.io as ragout_io

if __name__ == "__main__":
    full_description = camsa.full_description_template.format(
        names=camsa.CAMSA_AUTHORS,
        affiliations=camsa.AFFILIATIONS,
        dummy=" ",
        tool="Converting Ragout coords formatted blocks for further CAMSA processing.",
        information="For more information refer to {docs}".format(docs=camsa.CAMSA_DOCS_URL),
        contact=camsa.CONTACT)
    full_description = "=" * 80 + "\n" + full_description + "=" * 80 + "\n"
    parser = configargparse.ArgParser(description=full_description, formatter_class=configargparse.RawTextHelpFormatter,
                                      default_config_files=[os.path.join(camsa.root_dir, "utils", "ragout", "ragout_coords2camsa_seqi.ini"),
                                                            os.path.join(camsa.root_dir, "logging.ini")])

    parser.add_argument("--version", action="version", version=camsa.VERSION)
    parser.add_argument("-o", "--output-file", metavar="OUTPUT", dest="output", type=configargparse.FileType("wt"), default=sys.stdout,
                        help="A file to which the CAMSA readable fragments lengths would be written. Standard extension is \".camsa.lengths\".\nDEFAULT: stdout")
    parser.add_argument("--o-format", type=str, help="")
    parser.add_argument("ragout_coords", type=str, help="A path to ragout coords file")
    parser.add_argument("--genomes-order", type=str, default="",
                        help="A coma separated list of genome names, used to determine the output order for each homologous synteny block.\nDEFAULT: \"\" (i.e., sorted list of good (all - bad) genomes)")
    parser.add_argument("--ann-genomes", type=str,
                        default="A coma separated list of genome names, to be used when annotation all synteny blocks in the output.\nDEFAULT: \"\" (i.e., equals to the --genomes-order value)")
    parser.add_argument("--ann-delimiter", type=str, default=";")
    parser.add_argument("--filter-indels", action="store_true", dest="filter_indels", default=False,
                        help="A flag to indicate whether to filter out all synteny blocks, that are not present exactly once in each of the good (all-bad) genomes.")
    parser.add_argument("--no-filter-indels", action="store_false", dest="filter_indels", default=False,
                        help="A flag to indicate whether to filter out all synteny blocks, that are not present exactly once in each of the good (all-bad) genomes.")
    parser.add_argument("--filter-duplications", action="store_true", dest="filter_duplications", default=True,
                        help="A flag to indicate whether to filter out all synteny blocks, that are present more than once in at least on of the good (all-bad) genomes.")
    parser.add_argument("--no-filter-duplications", action="store_false", dest="filter_duplications", default=True,
                        help="A flag to indicate whether to filter out all synteny blocks, that are present more than once in at least on of the good (all-bad) genomes.")
    parser.add_argument("--good-genomes", type=str, default="", help="A coma separated list of genome names, to be processed and conversed.\nDEFAULT: \"\" (i.e., all genomes are good)")
    parser.add_argument("--bad-genomes", type=str, default="", help="A coma separated list of genome names, to be excluded from processing and conversion.\nDEFAULT: \"\" (i.e., no genomes are bad)")
    parser.add_argument("--c-logging-level", dest="logging_level", default=logging.INFO, type=int,
                        choices=[logging.NOTSET, logging.DEBUG, logging.INFO, logging.WARNING, logging.ERROR, logging.CRITICAL],
                        help="Logging level for the converter.\nDEFAULT: {info}".format(info=logging.INFO))
    parser.add_argument("--c-logging-formatter-entry",
                        help="Format string for python logger.")
    parser.add_argument("--o-delimiter", type=str, default="\t",
                        help="A single character string, used as a delimiter in the output (t)/(c)sv file.\nDEFAULT: \\t")
    args = parser.parse_args()
    start_time = datetime.datetime.now()

    logger = logging.getLogger("CAMSA.utils.ragout_coords2camsa_seqi")
    ch = logging.StreamHandler()
    ch.setLevel(args.logging_level)
    logger.setLevel(args.logging_level)
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
    if args.genomes_order != "":
        genomes_order = args.genomes_order.split(",")
    else:
        genomes_order = sorted(all_filtered_genomes)
    for genome in genomes_order:
        if genome not in all_filtered_genomes:
            logger.critical("Genome \"{genome}\" was specified with teh --genomes-order flag but was not found in the set of filtered genomes \{{filtered_genomes}\}"
                            "".format(genome=genome, filtered_genomes=",".join(sorted(all_filtered_genomes))))
            exit(1)
    if args.ann_genomes == "":
        args.ann_genomes = genomes_order
    else:
        args.ann_genomes = args.ann_genomes.split(",")
        for genome in args.ann_genomes:
            if genome not in genomes_order:
                logger.critical("Annotation genomes {ann_genome} was not present in reference genomes {ref_genomes}"
                                "".format(ann_genome=args.genome, ref_genomes=",".join(sorted(all_filtered_genomes))))
                exit(1)
    sequences = []
    for block_id in sorted(blocks_by_ids.keys()):
        blocks_by_id = blocks_by_ids[block_id]
        blocks = []
        for genome in genomes_order:
            blocks.extend([block for block in blocks_by_id if block.parent_seq.genome_name == genome])
        annotations = []
        for ann_genome in args.ann_genomes:
            for block in blocks:
                if block.parent_seq.genome_name == ann_genome:
                    annotations.append(block.annotation_name)
        annotation = args.ann_delimiter.join(annotations)
        for block in blocks:
            sequences.append(Sequence(name=block.name, parent_seq_id=block.parent_seq.seq_name, start=block.start, end=block.end, strand=block.strand, annotation=annotation, seq_group_id=block.parent_seq.genome_name))
    write_seqi(sequences=sequences, destination=args.output, output_setup=args.o_format, delimiter=args.o_delimiter)
    logger.info("All done!")
    end_time = datetime.datetime.now()
    logger.info("Elapsed time: {el_time}".format(el_time=str(end_time - start_time)))
