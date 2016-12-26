#! /usr/bin/env python
# -*- coding: utf-8 -*-
import datetime
import os
import sys
from collections import defaultdict

import configargparse
import logging
from bg.grimm import GRIMMReader

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))))

import camsa
from camsa.core.data_structures import AssemblyPoint
import camsa.core.io as camsa_io


def get_assembly_points(genomes):
    result = []
    for genome, chromosomes in genomes.items():
        for chr_type, blocks in chromosomes:
            left_blocks, right_blocks = blocks[:-1], blocks[1:]
            for l_block, r_block in zip(left_blocks, right_blocks):
                l_block_sign, l_block_name = l_block
                r_block_sign, r_block_name = r_block
                ap = AssemblyPoint(seq1=l_block_name, seq2=r_block_name, seq1_or=l_block_sign, seq2_or=r_block_sign, sources=[genome])
                result.append(ap)
            if chr_type == "@":
                first_block_sign, first_block_name = blocks[0]
                last_block_sign, last_block_name = blocks[-1]
                ap = AssemblyPoint(seq1=last_block_name, seq1_or=last_block_sign, seq2=first_block_name, seq2_or=first_block_sign, sources=[genome])
                result.append(ap)
    return result


def get_blocks_cnt(genomes):
    block_counter = defaultdict(lambda: defaultdict(int))
    for genome, chromosomes in genomes.items():
        for _, blocks in chromosomes:
            for block_sign, block_name in blocks:
                block_counter[block_name][genome] += 1
    return block_counter


def remove_bad_blocks(genomes, bad_blocks):
    for genome_name in list(genomes.keys()):
        new_chromosomes = []
        for chr_type, old_blocks in genomes[genome_name]:
            new_blocks = []
            for block_sign, block_name in old_blocks:
                if block_name in bad_blocks:
                    continue
                new_blocks.append((block_sign, block_name))
            if len(new_blocks) != 0:
                new_chromosomes.append((chr_type, new_blocks))
        genomes[genome_name] = new_chromosomes


def remove_duplications(genomes):
    blocks_cnt = get_blocks_cnt(genomes=genomes)
    bad_blocks = set()
    for block_name, block_cnt in blocks_cnt.items():
        cnt_per_genomes = block_cnt.values()
        if any(map(lambda x: x > 1, cnt_per_genomes)):
            bad_blocks.add(block_name)
    remove_bad_blocks(genomes=genomes, bad_blocks=bad_blocks)


def remove_indels(genomes):
    blocks_cnt = get_blocks_cnt(genomes=genomes)
    bad_blocks = set()
    all_genomes = set(genomes.keys())
    for block_name, block_cnt in blocks_cnt.items():
        blocks_genomes = set(block_cnt.keys())
        if len(all_genomes - blocks_genomes) != 0:
            bad_blocks.add(block_name)
    remove_bad_blocks(genomes=genomes, bad_blocks=bad_blocks)

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
    parser.add_argument("--trimmer-char", default=".", type=str, help="A character, which first occurrence in the genome name indicates a position for trimming.\nDEFAULT: .")
    parser.add_argument("--good-genomes", type=str, default="", help="A coma separated list of genome names, to be processed and conversed.\nDEFAULT: \"\" (i.e., all genomes are good)")
    parser.add_argument("--bad-genomes", type=str, default="", help="A coma separated list of genome names, to be excluded from processing and conversion.\nDEFAULT: \"\" (i.e., no genomes are bad)")
    parser.add_argument("--filter-duplications", dest="filter_duplications", default=False, action="store_true", help="")
    parser.add_argument("--filter-indels", dest="filter_indels", default=False, action="store_true", help="")
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

    genomes = defaultdict(list)
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
                        current_genome = current_genome.split(args.trimmer_char, 1)[0]
                elif current_genome is not None:
                    current_chromosome = []
                    chr_type, blocks = GRIMMReader.parse_data_string(data_string=line)
                    genomes[current_genome].append((chr_type, blocks))
    if args.good_genomes != "":
        good_genomes = args.good_genomes.split(",")
        if args.trim_names:
            good_genomes = [genome_name.split(args.trimmer_char, 1)[0] for genome_name in good_genomes]
        for genome_name in list(genomes.keys()):
            if genome_name not in good_genomes:
                del genomes[genome_name]
    if args.bad_genomes != "":
        bad_genomes = args.bad_genomes.split(",")
        if args.trim_names:
            bad_genomes = [genome_name.split(args.trimmer_char, 1)[0] for genome_name in bad_genomes]
        for genome_name in list(genomes.keys()):
            if genome_name in bad_genomes:
                del genomes[genome_name]
    if args.filter_duplications:
        remove_duplications(genomes=genomes)
    result = get_assembly_points(genomes=genomes)
    logger.info("Writing output to file \"{file_name}\"".format(file_name=args.output))
    camsa_io.write_assembly_points(assembly_points=result, destination=args.output, output_setup=args.o_format, delimiter=args.o_delimiter)
    logger.info("Elapsed time: {el_time}".format(el_time=str(datetime.datetime.now() - start_time)))






