#! /usr/bin/env python
# -*- coding: utf-8 -*-
import csv
import datetime
import logging
import os
import sys

import configargparse
from Bio import SeqIO


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
    parser.add_argument("--ref-genomes", type=str, default="")
    parser.add_argument("--ann-genomes", type=str, default="")
    parser.add_argument("--ann-delimiter", type=str, default=";")
    parser.add_argument("--filter-indels", action="store_true", dest="filter_indels", default=False)
    parser.add_argument("--filter-duplications", action="store_true", dest="filter_duplications", default=False)
    parser.add_argument("--good-genomes", type=str, default="", help="A coma separated list of genome names, to be processed and conversed.\nDEFAULT: \"\" (i.e., all genomes are good)")
    parser.add_argument("--bad-genomes", type=str, default="", help="A coma separated list of genome names, to be excluded from processing and conversion.\nDEFAULT: \"\" (i.e., no genomes are bad)")
    parser.add_argument("--sbs", action="store_true", default=False, dest="silent_block_skip")
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

    if args.ref_genomes == "":
        logger.info("Reference genomes were not specified")
        args.ref_genomes = sorted(all_filtered_genomes)
        logger.info("Setting reference genomes to {ref_genome}".format(ref_genome=args.ref_genome))
    else:
        args.ref_genomes = args.ref_genomes.split(",")
        for genome in args.ref_genomes:
            if genome not in all_filtered_genomes:
                logger.critical("Reference genome {ref_genome} was not present in all filtered genome {filtered_genomes}"
                                "".format(ref_genome=args.genome, filtered_genomes=",".join(all_filtered_genomes)))
                exit(1)
    if args.ann_genomes == "":
        args.ann_genomes = args.ref_genomes
    else:
        args.ann_genomes = args.ann_genomes.split(",")
        for genome in args.ann_genomes:
            if genome not in args.ref_genomes:
                logger.critical("Annotation genomes {ann_genome} was not present in reference genomes {ref_genomes}"
                                "".format(ann_genome=args.genome, ref_genomes=",".join(all_filtered_genomes)))
                exit(1)
    sequences = []
    for block_id in sorted(blocks_by_ids.keys()):
        blocks = blocks_by_ids[block_id]
        blocks_by_ref_genomes = [block for block in blocks if block.parent_seq.genome_name in args.ref_genomes]
        if len(blocks_by_ref_genomes) == 0:
            if not args.silent_block_skip:
                logger.critical("For blocks with id {block_id} not a single instance was present in the reference genome {ref_genome}. Flag \"--sbs\" was not set, and this event is thus critical."
                                "".format(block_id=block_id, ref_genome=args.ref_genome))
            else:
                logger.warning("For blocks with id {block_id} not a single instance was present in the reference genome \"{ref_genome}\". Flag \"--sbs\" was set, thus silently ignoring this case."
                               "".format(block_id=block_id, ref_genome=args.ref_genome))
                continue
        # if len(blocks_by_ref_genomes) > 1:
        #     logger.warning("More than a single block with id {block_id} was found in the reference genome \"{ref_genome}\". Randomly choosing one such block (shall not be a problem, as they must be merely identical)"
        #                    "".format(block_id=block_id, ref_genome=args.ref_genome))
        blocks = []
        for ref_genome in args.ref_genomes:
            for block in blocks_by_ref_genomes:
                if block.parent_seq.genome_name == ref_genome:
                    blocks.append(block)
        annotations = []
        for ann_genome in args.ann_genomes:
            for block in blocks_by_ref_genomes:
                if block.parent_seq.genome_name == ann_genome:
                    annotations.append(block.annotation_name)
        annotation = args.ann_delimiter.join(annotations)
        block = blocks[0]
        sequences.append(Sequence(name=block.name, parent_seq_id=block.parent_seq.seq_name, start=block.start, end=block.end, strand=block.strand, annotation=annotation))

    write_seqi(sequences=sequences, destination=args.output, output_setup=args.o_format, delimiter=args.o_delimiter)
    logger.info("All done!")
    end_time = datetime.datetime.now()
    logger.info("Elapsed time: {el_time}".format(el_time=str(end_time - start_time)))