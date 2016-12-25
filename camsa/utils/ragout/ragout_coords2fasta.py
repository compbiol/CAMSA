# -*- coding: utf-8 -*-
import datetime
import logging
import os
import sys
from collections import defaultdict

import configargparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from utils.ragout.shared import filter_blocks_by_good_genomes, filter_blocks_by_bad_genomes, get_all_genomes_from_blocks

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))))

import camsa
import camsa.utils.ragout.io as ragout_io
from camsa.utils.ragout.shared import filter_indels, filter_duplications

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
                                      default_config_files=[os.path.join(camsa.root_dir, "logging.ini")])
    parser.add_argument("-c", "--config", is_config_file=True, help="Config file overwriting some of the default settings as well as any flag starting with \"--\".")

    parser.add_argument("--version", action="version", version=camsa.VERSION)
    parser.add_argument("ragout_coords", type=str, help="A path to ragout coords file")
    parser.add_argument("--ref-genome", type=str, default="")
    parser.add_argument("--filter-indels", type=bool, action="store_true", dest="filter_indels", default=False)
    parser.add_argument("--filter_duplications", type=bool, action="store_true", dest="filter_duplications", default=False)
    parser.add_argument("--good-genomes", type=str, default="", help="A coma separated list of genome names, to be processed and conversed.\nDEFAULT: \"\" (i.e., all genomes are good)")
    parser.add_argument("--bad-genomes", type=str, default="", help="A coma separated list of genome names, to be excluded from processing and conversion.\nDEFAULT: \"\" (i.e., no genomes are bad)")
    parser.add_argument("--sbs", action="store_true", default=False, dest="silent_block_skip", type=bool)
    parser.add_argument("fasta", nargs="+")
    parser.add_argument("-o", "--output", type=configargparse.FileType("wt"), default=sys.stdout)
    parser.add_argument("--c-logging-level", dest="logging_level", default=logging.INFO, type=int,
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
        good_genomes = set(args.good_genomes.split(","))
        filter_blocks_by_good_genomes(blocks_by_ids=blocks_by_ids, good_genomes=good_genomes)
    if args.bad_genomes != "":
        bad_genomes = set(args.bad_genomes.split(","))
        filter_blocks_by_bad_genomes(blocks_by_ids=blocks_by_ids, bad_genomes=bad_genomes)
    if args.filter_indels:
        filter_indels(blocks_by_ids=blocks_by_ids, all_genomes_as_set=all_genomes)
    if args.filter_duplications:
        filter_duplications(blocks_by_ids=blocks_by_ids)
    all_filtered_genomes = get_all_genomes_from_blocks(blocks_as_ids=blocks_by_ids)
    if args.ref_genome == "":
        logger.info("Reference genome was not specified")
        args.ref_genome = sorted(all_filtered_genomes)[0]
        logger.info("Setting reference genome to {ref_genome}".format(ref_genome=args.ref_genome))
    else:
        if args.ref_genome not in all_filtered_genomes:
            logger.critical("Reference genome {ref_genome} was not present in all filtered genome {filtered_genomes}"
                            "".format(ref_genome=args.ref_genome, filtered_genomes=",".join(all_filtered_genomes)))
            exit(1)
    blocks_to_convert = []
    for block_id in sorted(blocks_by_ids.keys()):
        blocks = blocks_by_ids[block_id]
        blocks_by_ref_genome = [block for block in blocks if block.parent_seq.genome_name == args.ref_genome]
        if len(blocks_by_ref_genome) == 0:
            if not args.silent_block_skip:
                logger.critical("For blocks with id {block_id} not a single instance was present in the reference genome {ref_genome}. Flag \"--sbs\" was not set, and this event is thus critical."
                                "".format(block_id=block_id, ref_genome=args.ref_genome))
            else:
                logger.warning("For blocks with id {block_id} not a single instance was present in the reference genome \"{ref_genome}\". Flag \"--sbs\" was set, thus silently ignoring this case."
                               "".format(block_id=block_id, ref_genome=args.ref_genome))
        if len(blocks_by_ref_genome) > 1:
            logger.warning("More than a single block with id {block_id} was found in the reference genome \"{ref_genome}\". Randomly choosing one such block (shall not be a problem, as they must be merely identical)"
                           "".format(block_id=block_id, ref_genome=args.ref_genome))
        block = blocks_by_ref_genome[0]
        blocks_to_convert.append(block)

    blocks_by_seq_ids = defaultdict(list)
    for block in blocks_to_convert:
        blocks_by_seq_ids[block.parent_seq.seq_name].append(block)

    for f in args.fasta:
        logger.info("Processing fasta file: \"{file_name}\"".format(file_name=f))
        cnt = 0
        for record in SeqIO.parse(f, "fasta"):
            seq_id = record.id
            if seq_id not in blocks_by_ids:
                continue
            current_blocks = blocks_by_ids[seq_id]
            for block in current_blocks:
                if block.strand == "+":
                    out_seq = record.seq[block.start: block.end]
                else:
                    out_seq = record.seq[block.start: block.end].reverse_complement()
                out_record = SeqRecord(seq=out_seq, id=block.name, description="")
                SeqIO.write(sequences=out_record, handle=args.output, format="fasta")

    logger.info("All done!")
    end_time = datetime.datetime.now()
    logger.info("Elapsed time: {el_time}".format(el_time=str(end_time - start_time)))
