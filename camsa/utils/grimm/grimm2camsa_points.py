# -*- coding: utf-8 -*-
import datetime
import os
import sys

import configargparse
import logging
from bg.grimm import GRIMMReader

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))))

import camsa
from camsa.core.data_structures import AssemblyPoint
import camsa.core.io as camsa_io

if __name__ == "__main__":
    full_description = camsa.full_description_template.format(
        names=camsa.CAMSA_AUTHORS,
        affiliations=camsa.AFFILIATIONS,
        dummy=" ",
        tool="Converting GRIMM formatted scaffolding results for further CAMSA processing.",
        information="For more information refer to {docs}".format(docs=camsa.CAMSA_DOCS_URL),
        contact=camsa.CONTACT)
    full_description = "=" * 80 + "\n" + full_description + "=" * 80 + "\n"
    parser = configargparse.ArgParser(description=full_description, formatter_class=configargparse.RawTextHelpFormatter,
                                      default_config_files=[os.path.join(camsa.root_dir, "utils", "grimm", "grimm2camsa_points.ini"),
                                                            os.path.join(camsa.root_dir, "logging.ini")])
    parser.add_argument("grimm", nargs="+", help="A list of paths to files, that in GRIMM format contain information about scaffold assemblies.")
    parser.add_argument("--o-format", type=str, help="")
    parser.add_argument("-o", "--output", type=configargparse.FileType("wt"), default=sys.stdout, help="Path to the file, where converted assembly points will be stored.\nDEFAULT: stdout")
    parser.add_argument("--o-delimiter", default="\t", type=str, help="")
    parser.add_argument("--no-trim-names", dest="trim_names", default=True, action="store_false",
                        help="A flag to indicate, that genome names from grimm file need not to be trimmed by the first \".\"\nDEFAULT: true")
    parser.add_argument("--c-logging-level", type=int,
                        choices=[logging.NOTSET, logging.DEBUG, logging.INFO, logging.WARNING, logging.ERROR, logging.CRITICAL],
                        help="Logging level for the converter.\nDEFAULT: {info}".format(info=logging.INFO))
    parser.add_argument("--c-logging-formatter-entry",
                        help="Format string for python logger.")

    args = parser.parse_args()

    start_time = datetime.datetime.now()

    logger = logging.getLogger("CAMSA.utils.grimm2camsa_points")
    ch = logging.StreamHandler()
    ch.setLevel(args.c_logging_level)
    logger.setLevel(args.c_logging_level)
    logger.addHandler(ch)
    logger.info(full_description)
    logger.info(parser.format_values())
    ch.setFormatter(logging.Formatter(args.c_logging_formatter_entry))
    logger.info("Starting the converting process")

    result = []
    for file_name in args.grimm:
        logger.info("Processing file \"{file_name}\"".format(file_name=file_name))
        with open(file_name, "rt") as source:
            current_genome = None
            for line in source:
                line = line.strip()
                if len(line) == 0 or GRIMMReader.is_comment_string(data_string=line):
                    continue
                if GRIMMReader.is_genome_declaration_string(data_string=line):
                    current_genome = GRIMMReader.parse_genome_declaration_string(data_string=line).name
                    if args.trim_names:
                        current_genome = current_genome.split(".", 1)[0]
                elif current_genome is not None:
                    chr_type, blocks = GRIMMReader.parse_data_string(data_string=line)
                    left_blocks, right_blocks = blocks[:-1], blocks[1:]
                    for l_block, r_block in zip(left_blocks, right_blocks):
                        l_block_sign, l_block_name = l_block
                        r_block_sign, r_block_name = r_block
                        ap = AssemblyPoint(seq1=l_block_name, seq2=r_block_name, seq1_or=l_block_sign, seq2_or=r_block_sign, sources=[current_genome])
                        result.append(ap)
    logger.info("Writing output to file \"{file_name}\"".format(file_name=args.output))
    camsa_io.write_assembly_points(assembly_points=result, destination=args.output, output_setup=args.o_format, delimiter=args.o_delimiter)
    logger.info("Elapsed time: {el_time}".format(el_time=str(datetime.datetime.now() - start_time)))






