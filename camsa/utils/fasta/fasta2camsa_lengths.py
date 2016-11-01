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

if __name__ == "__main__":
    full_description = camsa.full_description_template.format(
        names=camsa.CAMSA_AUTHORS,
        affiliations=camsa.AFFILIATIONS,
        dummy=" ",
        tool="Converting FASTA formatted scaffolding results for further CAMSA processing.",
        information="For more information refer to {docs}".format(docs=camsa.CAMSA_DOCS_URL),
        contact=camsa.CONTACT)
    full_description = "=" * 80 + "\n" + full_description + "=" * 80 + "\n"
    parser = configargparse.ArgParser(description=full_description, formatter_class=configargparse.RawTextHelpFormatter)
    parser.add_argument("--version", action="version", version=camsa.VERSION)
    parser.add_argument("-o", "--output-file", metavar="OUTPUT", dest="output", type=configargparse.FileType("wt"), default=sys.stdout,
                        help="A file to which the CAMSA readable fragments lengths would be written. Standard extension is \".camsa.lengths\".\nDEFAULT: stdout")
    parser.add_argument("contigs", nargs="+", metavar="CONTIGS", type=configargparse.FileType("rt"), default=sys.stdin,
                        help="A list of input *.fasta files with contigs.\nDEFAULT: stdin")
    parser.add_argument("--c-logging-level", dest="logging_level", default=logging.INFO, type=int,
                        choices=[logging.NOTSET, logging.DEBUG, logging.INFO, logging.WARNING, logging.ERROR, logging.CRITICAL],
                        help="Logging level for the converter.\nDEFAULT: {info}".format(info=logging.INFO))
    parser.add_argument("--c-logging-formatter-entry",
                        help="Format string for python logger.")
    args = parser.parse_args()
    start_time = datetime.datetime.now()

    logger = logging.getLogger("CAMSA.utils.fasta2camsa_lengths")
    ch = logging.StreamHandler()
    ch.setLevel(args.logging_level)
    logger.setLevel(args.logging_level)
    logger.addHandler(ch)
    logger.info(full_description)
    logger.info(parser.format_values())
    ch.setFormatter(logging.Formatter(args.c_logging_formatter_entry))
    logger.info("Starting the converting process")

    entries = {}
    for f in args.contigs:
        logger.info("Processing file: \"{file_name}\"".format(file_name=f))
        cnt = 0
        for record in SeqIO.parse(f, "fasta"):
            entries[record.id] = len(record.seq)
            cnt += 1
        logger.info("Processed {cnt} fasta records".format(cnt=cnt))
    writer = csv.writer(args.output)
    logger.info("Writing the CAMSA lengths output")
    writer.writerow(["seq_id", "seq_length"])
    for key in sorted(entries.keys()):
        writer.writerow([key, entries[key]])
    logger.info("All done!")
    end_time = datetime.datetime.now()
    logger.info("Elapsed time: {el_time}".format(el_time=str(end_time - start_time)))
