#! /usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
import csv
import datetime
import logging
import sys

from Bio import SeqIO

if __name__ == "__main__":
    full_description = "=" * 80 + \
                       "\nSergey Aganezov & Max A. Alekseyev (c)\n" + \
                       "Computational Biology Institute, The George Washington University.\n\n" + \
                       "Extracting contigs lengths from *.fasta files for further CAMSA processing.\n" + \
                       "With any questions, please, contact Sergey Aganezov [aganezov(at)gwu.edu].\n" + \
                       "=" * 80 + "\n"
    parser = argparse.ArgumentParser(description=full_description, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-o", "--output-file", meta="OUTPUT", dest="output", type=argparse.FileType("wt"), default=sys.stdout,
                        help="A file to which the CAMSA readable fragments lengths would be written. Standard extension is \".camsa.lengths\".\nDEFAULT: stdout")
    parser.add_argument("contigs", nargs="+", meta="CONTIGS", dest="contig_files", type=argparse.FileType("rt"), default=sys.stdin,
                        help="A list of input *.fasta files with contigs.\nDEFAULT: stdin")
    parser.add_argument("--logging-level", dest="logging_level", default=logging.INFO, type=int,
                        choices=[logging.NOTSET, logging.DEBUG, logging.INFO, logging.WARNING, logging.ERROR, logging.CRITICAL],
                        help="Logging level for the converter.\nDEFAULT: {info}".format(info=logging.INFO))
    args = parser.parse_args()
    start_time = datetime.datetime.now()

    logger = logging.getLogger("fasta2camsa_lengths")
    ch = logging.StreamHandler()
    ch.setLevel(args.logging_level)
    logger.setLevel(args.logging_level)
    logger.addHandler(ch)
    logger.info(full_description)
    logger.info("Starting the converting process")

    entries = {}
    for f in args.contig_files:
        logger.info("Processing file: \"{file_name}\"".format(file_name=f))
        cnt = 0
        for record in SeqIO.parse(f, "fasta"):
            entries[record.id] = len(record.seq)
            cnt += 1
        logger.info("Processed {cnt} fasta records".format(cnt=cnt))
    writer = csv.writer(args.output)
    logger.info("Writing the CAMSA lengths output")
    writer.writerow(["ctg_id", "ctg_length"])
    for key in sorted(entries.keys()):
        writer.writerow([key, entries[key]])
    logger.info("All done!")
    end_time = datetime.datetime.now()
    logger.info("Elapsed time: {el_time}".format(el_time=str(end_time - start_time)))
