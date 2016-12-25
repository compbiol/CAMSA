# -*- coding: utf-8 -*-
import os
import sys

import configargparse
import logging

from utils.ragout.shared import filter_blocks_by_good_genomes, filter_blocks_by_bad_genomes

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
    # parser.add_argument("fasta", nargs="+")
    parser.add_argument("--c-logging-level", dest="logging_level", default=logging.INFO, type=int,
                        choices=[logging.NOTSET, logging.DEBUG, logging.INFO, logging.WARNING, logging.ERROR, logging.CRITICAL],
                        help="Logging level for the converter.\nDEFAULT: {info}".format(info=logging.INFO))
    parser.add_argument("--c-logging-formatter-entry",
                        help="Format string for python logger.")
    args = parser.parse_args()

    sequences_by_ids, blocks_by_ids = ragout_io.read_from_file(path=args.ragout_coords, silent_fail=False, delimiter="\t")
    if args.good_genomes != "":
        good_genomes = set(args.good_genomes.split(","))
        filter_blocks_by_good_genomes(blocks_by_ids=blocks_by_ids, good_genomes=good_genomes)
    if args.bad_genomes != "":
        bad_genomes = set(args.bad_genomes.split(","))
        filter_blocks_by_bad_genomes(blocks_by_ids=blocks_by_ids, bad_genomes=bad_genomes)
    if args.filter_indels:
        filter_indels(blocks_by_ids=blocks_by_ids)
    if args.filter_duplications:
        filter_duplications(blocks_by_ids=blocks_by_ids)

    a = 5
