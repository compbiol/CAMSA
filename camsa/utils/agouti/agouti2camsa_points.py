#! /usr/bin/env python
# -*- coding: utf-8 -*-

import datetime
import os
import sys
import more_itertools
from collections import defaultdict

import configargparse
import logging

import camsa
from camsa.core.data_structures import AssemblyPoint
import camsa.core.io as camsa_io

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))))


def get_assembly_points(agouti_path, source):
    result = []
    for left, right in more_itertools.windowed(agouti_path, n=2):
        ap = AssemblyPoint(seq1=left, seq2=right, seq1_or="?", seq2_or="?", sources=[str(source)])
        result.append(ap)
    return result


if __name__ == "__main__":
    full_description = camsa.full_description_template.format(
        names=camsa.CAMSA_AUTHORS,
        affiliations=camsa.AFFILIATIONS,
        dummy=" ",
        tool="Converting AGOUTI formatted scaffolding results for further CAMSA processing.",
        information="For more information refer to {docs}".format(docs=camsa.CAMSA_DOCS_URL),
        contact=camsa.CONTACT)
    full_description = "=" * 80 + "\n" + full_description + "=" * 80 + "\n"
    parser = configargparse.ArgParser(description=full_description, formatter_class=configargparse.RawTextHelpFormatter,
                                      default_config_files=[os.path.join(camsa.root_dir, "utils", "agouti", "agouti2camsa_points.ini"),
                                                            os.path.join(camsa.root_dir, "logging.ini")])
    parser.add_argument("agouti", nargs="+", help="A list of paths to files, that in AGOUTI format contain information about scaffold assemblies.")
    parser.add_argument("--o-format", type=str, help="")
    parser.add_argument("-o", "--output", type=configargparse.FileType("wt"), default=sys.stdout, help="Path to the file, where converted assembly points will be stored.\nDEFAULT: stdout")
    parser.add_argument("--o-delimiter", default="\t", type=str, help="")
    parser.add_argument("--source", default=None, help="A value to be used in the \"source\" column in the output.\n"
                                                       " If not specified, the name of the file for each set of paths will be used for all assembly points inferred from the corresponding paths.")
    parser.add_argument("--i-delimiter", default=",", type=str, help="A character to be used as a delimiter during the parsing of the AGOUTI data sets")
    parser.add_argument("--c-logging-level", type=int,
                        choices=[logging.NOTSET, logging.DEBUG, logging.INFO, logging.WARNING, logging.ERROR, logging.CRITICAL],
                        help="Logging level for the converter.\nDEFAULT: {info}".format(info=logging.INFO))
    parser.add_argument("--c-logging-formatter-entry",
                        help="Format string for python logger.")
    args = parser.parse_args()

    start_time = datetime.datetime.now()

    logger = logging.getLogger("CAMSA.utils.agouti2camsa_points")
    ch = logging.StreamHandler()
    ch.setLevel(args.c_logging_level)
    logger.setLevel(args.c_logging_level)
    logger.addHandler(ch)
    logger.info(full_description)
    logger.info(parser.format_values())
    ch.setFormatter(logging.Formatter(args.c_logging_formatter_entry))
    logger.info("Starting the converting process")

    paths = []
    for file_name in args.agouti:
        logger.info("Processing file \"{file_name}\"".format(file_name=file_name))
        with open(file_name, "rt") as source:
            for line in source:
                line = line.strip()
                if len(line) == 0 or line.startswith("#") or line.startswith(">"):
                    continue
                path = line.split(args.i_delimiter)
                if args.source is not None:
                    source = args.source
                else:
                    source = os.path.basename(file_name)
                    source = os.path.splitext(source)[0]
                paths.append((source, path))
    logger.info("A total of {paths_cnt} paths were extracted from {file_cnt} files"
                "".format(paths_cnt=len(paths), file_cnt=len(args.agouti)))
    assembly_points = []
    for source, path in paths:
        if len(path) <= 1:
            logger.warning("Encountered a path of length <= 1 {{{path}}}; skipping"
                           "".format(path=",".join(path)))
            continue
        assembly_points.extend(get_assembly_points(agouti_path=path, source=source))
    logger.info("Writing output to file \"{file_name}\"".format(file_name=args.output))
    camsa_io.write_assembly_points(assembly_points=assembly_points, destination=args.output, output_setup=args.o_format, delimiter=args.o_delimiter)
    logger.info("Elapsed time: {el_time}".format(el_time=str(datetime.datetime.now() - start_time)))
