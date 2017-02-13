#! /usr/bin/env python
# -*- coding: utf-8 -*-
import datetime
import logging
import os
import sys
from collections import defaultdict

import configargparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

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
        tool="Converting Ragout formatted blocks into FASTA format.",
        information="For more information refer to {docs}".format(docs=camsa.CAMSA_DOCS_URL),
        contact=camsa.CONTACT)
    full_description = "=" * 80 + "\n" + full_description + "=" * 80 + "\n"
    parser = configargparse.ArgParser(description=full_description, formatter_class=configargparse.RawTextHelpFormatter,
                                      default_config_files=[os.path.join(camsa.root_dir, "logging.ini")])
    parser.add_argument("-c", "--config", is_config_file=True, help="Config file overwriting some of the default settings as well as any flag starting with \"--\".")

    parser.add_argument("--version", action="version", version=camsa.VERSION)
    parser.add_argument("ragout_coords", type=str, help="A path to ragout coords file")
    parser.add_argument("--ann-genomes", type=str, default="")
    parser.add_argument("--ann-delimiter", type=str, default=";")
    parser.add_argument("--filter-indels", action="store_true", dest="filter_indels", default=False)
    parser.add_argument("--filter-duplications", action="store_true", dest="filter_duplications", default=False)
    parser.add_argument("--good-genomes", type=str, default="", help="A coma separated list of genome names, to be processed and conversed.\nDEFAULT: \"\" (i.e., all genomes are good)")
    parser.add_argument("--bad-genomes", type=str, default="", help="A coma separated list of genome names, to be excluded from processing and conversion.\nDEFAULT: \"\" (i.e., no genomes are bad)")
    parser.add_argument("fasta", nargs="+")
    parser.add_argument("-o", "--output", type=configargparse.FileType("wt"), default=sys.stdout)
    parser.add_argument("--o-genomes", dest="ref_genomes", type=str, default="a string of coma separated names of genomes, who will")
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
    if args.ref_genomes == "":
        logger.info("Output genomes were not specified")
        args.ref_genomes = sorted(all_filtered_genomes)
        logger.info("Setting output genomes to {ref_genome}".format(ref_genome=args.ref_genomes))
    else:
        args.ref_genomes = args.ref_genomes.split(",")
        for genome in args.ref_genomes:
            if genome not in all_filtered_genomes:
                logger.critical("Output genome {ref_genome} was not present in all filtered genomes {filtered_genomes}"
                                "".format(ref_genome=args.genome, filtered_genomes=",".join(all_filtered_genomes)))
                exit(1)
    if args.ann_genomes == "":
        args.ann_genomes = args.ref_genomes
    else:
        args.ann_genomes = args.ann_genomes.split(",")
        for genome in args.ann_genomes:
            if genome not in all_filtered_genomes:
                logger.critical("Annotation genomes {ann_genome} was not present in all filtered genomes {filtered_genomes}"
                                "".format(ann_genome=args.genome, filtered_genomes=",".join(all_filtered_genomes)))
                exit(1)
    blocks_to_convert = []
    for block_id in sorted(blocks_by_ids.keys()):
        blocks = blocks_by_ids[block_id]
        blocks_by_ref_genomes = [block for block in blocks if block.parent_seq.genome_name in args.ref_genomes]
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
        block = blocks[0]
        if len(annotations) != 0:
            block._annotation = args.ann_delimiter.join(annotations)
        blocks_to_convert.append(block)

    blocks_by_seq_ids = defaultdict(list)
    for block in blocks_to_convert:
        blocks_by_seq_ids[block.parent_seq.seq_name].append(block)

    processed = set()
    for f in args.fasta:
        logger.info("Processing fasta file: \"{file_name}\"".format(file_name=f))
        cnt = 0
        for record in SeqIO.parse(f, "fasta"):
            seq_id = record.id
            if seq_id not in blocks_by_seq_ids:
                continue
            current_blocks = blocks_by_seq_ids[seq_id]
            current_blocks = [block for block in current_blocks if block.name not in processed]
            for block in current_blocks:
                if block.strand == "+":
                    out_seq = record.seq[block.start: block.end]
                else:
                    out_seq = record.seq[block.start: block.end].reverse_complement()
                out_record = SeqRecord(seq=out_seq, id=str(block.name), description=block.annotation_name)
                SeqIO.write(sequences=out_record, handle=args.output, format="fasta")
                processed.add(block.name)

    logger.info("All done!")
    end_time = datetime.datetime.now()
    logger.info("Elapsed time: {el_time}".format(el_time=str(end_time - start_time)))
