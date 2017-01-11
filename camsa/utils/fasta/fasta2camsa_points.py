#! /usr/bin/env python
# -*- coding: utf-8 -*-
import csv
import datetime
import logging
import os
import subprocess
import sys
from collections import defaultdict, deque
from os.path import isfile

import configargparse

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))))

import camsa


class CoordsEntry(object):
    __slots__ = ['scaffold_name', 'fragment_name', 'scaffold_start', 'scaffold_end',
                 'fragment_start', 'fragment_end', 'fragment_length', 'fragment_coverage']

    def __init__(self, scaffold_name, fragment_name, scaf_start, scaf_end, frag_start, frag_end, frag_length, frag_cov):
        self.scaffold_name = scaffold_name
        self.fragment_name = fragment_name
        self.scaffold_start = scaf_start
        self.scaffold_end = scaf_end
        self.fragment_start = frag_start
        self.fragment_end = frag_end
        self.fragment_length = frag_length
        self.fragment_coverage = frag_cov

    @property
    def fragment_orientation(self):
        return "+" if self.fragment_start < self.fragment_end else "-"

    @property
    def mid_coordinate(self):
        return (self.scaffold_start + self.scaffold_end) / 2


class CoordsParser(object):
    """
    This is a simple parser, that always anticipates, that entries coordinates would be in two first meta-columns

    entries names would be in the last
    length and coverage of the query entry is anticipated to be in the second and third columns from the end, respectively
    """

    def __init__(self, source, skip_lines=3, has_header=True, lower_frag_cover=90.0, collapse_consecutive=True):
        self.source = source
        self.skip_lines = skip_lines
        self.has_header = has_header
        self.lower_frag_cover = lower_frag_cover
        self.collapse_consecutive = collapse_consecutive

    def parse_data(self):
        result = defaultdict(list)
        for cnt, entry in enumerate(self.source):
            if cnt < self.skip_lines:
                continue
            elif cnt in [self.skip_lines, self.skip_lines + 1] and self.has_header:
                continue
            coords_entry = self.from_string(string=entry)
            # if coords_entry.fragment_coverage < self.lower_frag_cover:
            #     continue
            result[coords_entry.scaffold_name].append(coords_entry)
        for name in result.keys():
            result[name] = sorted(result[name], key=lambda coord_entry: coords_entry.mid_coordinate)
        if self.collapse_consecutive:
            for name in list(result.keys()):
                entries = result[name]
                new_entries = []
                current_entry = entries[0]
                for entry in entries[1:]:
                    if entry.fragment_name == current_entry.fragment_name \
                            and entry.fragment_orientation == current_entry.fragment_orientation:
                        current_entry.fragment_coverage += entry.fragment_coverage
                        if entry.fragment_orientation == "+":
                            current_entry.scaffold_end = entry.scaffold_end
                        else:
                            current_entry.scaffold_start = entry.scaffold_start
                    else:
                        new_entries.append(current_entry)
                        current_entry = entry
                new_entries.append(current_entry)
                result[name] = new_entries
        for name in list(result.keys()):
            entries = result[name]
            new_entries = []
            for entry in entries:
                if entry.fragment_coverage >= self.lower_frag_cover:
                    new_entries.append(entry)
            result[name] = new_entries
        return result

    @staticmethod
    def from_string(string):
        data = [entry.strip() for entry in string.split(" | ")]
        reference_coords_data, query_coords_data = data[0], data[1]
        names_data = data[-1]
        lengths_data, coverage_data = data[-3], data[-2]
        ref_start, ref_end = [int(entry.strip()) for entry in reference_coords_data.split()]
        frag_start, frag_end = [int(entry.strip()) for entry in query_coords_data.split()]
        ref_name, frag_name = [entry.strip() for entry in names_data.split()]
        frag_length = float(lengths_data.split()[1].strip())
        frag_cov = float(coverage_data.split()[1].strip())
        return CoordsEntry(scaffold_name=ref_name, fragment_name=frag_name,
                           scaf_start=ref_start, scaf_end=ref_end,
                           frag_start=frag_start, frag_end=frag_end,
                           frag_length=frag_length, frag_cov=frag_cov)


def get_file_prefix(file_name):
    return os.path.basename(os.path.splitext(file_name)[0])


def run_nucmer(contigs_file_name, reference_file_name, nucmer_executable_path, cli_arguments, output_dir, logs_dir, nucmer_stdout=None, nucmer_stderr=None):
    prefix = get_file_prefix(reference_file_name)
    if nucmer_stdout is None:
        nucmer_stdout = os.path.join(logs_dir, "nucmer_{prefix}.stdout.txt".format(prefix=prefix))
    if nucmer_stderr is None:
        nucmer_stderr = os.path.join(logs_dir, "nucmer_{prefix}.stderr.txt".format(prefix=prefix))
    cli_args = [nucmer_executable_path] + cli_arguments.split(" ") + ['-p', os.path.join(output_dir, prefix), reference_file_name, contigs_file_name, '>', nucmer_stdout, '2>', nucmer_stderr]
    logger.info("\t" + " ".join(cli_args))
    exitcode = subprocess.call(" ".join(cli_args), shell=True)
    if exitcode == 0:
        logger.info("NUCmer finished running for \"{scaffolds_file}\" scaffolds file.")
        logger.debug("NUCmer alignment output is stored in \"{delta_output}\" file."
                     "".format(delta_output=os.path.join(output_dir, "{prefix}.delta".format(prefix=prefix))))
        logger.debug("NUCmer log is stored in \"{nucmer_log}\".".format(scaffolds_file=reference_file_name, nucmer_log=nucmer_stdout))
    else:
        logger.error("NUCmer exited with non-zero code, running for \"{scaffolds}\" scaffolds file.".format(scaffolds=reference_file_name))
        logger.error("NUCmer logs are stored in:")
        logger.error("\tstdout: \"{nucmer_log}\"".format(scaffolds=reference_file_name, nucmer_log=nucmer_stdout, nucmer_err=nucmer_stderr))
        logger.error("\tstderr: \"{nucmer_err}\"".format(scaffolds=reference_file_name, nucmer_log=nucmer_stdout, nucmer_err=nucmer_stderr))
    return exitcode


def run_show_coords(delta_file_name, output_dir, logs_dir, show_coords_executable_path, cli_arguments, show_coords_stderr=None):
    prefix = get_file_prefix(delta_file_name)
    if prefix.endswith(".filtered"):
        prefix = prefix[:-9]
    output_file_name = os.path.join(output_dir, prefix + ".coords")
    if show_coords_stderr is None:
        show_coords_stderr = os.path.join(logs_dir, "show_coords_{prefix}.stderr.txt".format(prefix=prefix))
    cli_args = [show_coords_executable_path] + cli_arguments.split(" ") + [delta_file_name, '>', output_file_name, '2>', show_coords_stderr]
    logger.info("\t" + " ".join(cli_args))
    exitcode = subprocess.call(" ".join(cli_args), shell=True)
    if exitcode == 0:
        logger.info("show-coords util has finished running for \"{delta_file}\".".format(delta_file=delta_file_name))
        logger.debug("show-coords output is stored in \"{coords_output}\" file.".format(coords_output=output_file_name))
    else:
        logger.error("show-coords exited with non-zero code, running for \"{delta_file}\".".format(delta_file=delta_file_name, show_coords_error_log=show_coords_stderr))
        logger.error("\tshow-coords error log is stored in \"{show_coords_error_log}\"".format(delta_file=delta_file_name, show_coords_error_log=show_coords_stderr))
    return exitcode


def run_delta_filter(delta_file_name, output_dir, logs_dir, delta_filter_executable_path, cli_arguments, delta_filter_stderr=None):
    prefix = get_file_prefix(delta_file_name)
    output_file_name = os.path.join(output_dir, prefix + ".filtered.delta")
    if delta_filter_stderr is None:
        delta_filter_stderr = os.path.join(logs_dir, "delta_filter_{prefix}.stderr.txt".format(prefix=prefix))
    cli_args = [delta_filter_executable_path] + cli_arguments.split(" ") + [delta_file_name, ">", output_file_name, "2>", delta_filter_stderr]
    logger.info("\t" + " ".join(cli_args))
    exitcode = subprocess.call(" ".join(cli_args), shell=True)
    if exitcode == 0:
        logger.info("delta-filter util has finished running for \"{delta_file}\"".format(delta_file=delta_file_name))
        logger.debug("delta-filter output is stored in \"{delta_filter_output}\" file".format(delta_filter_output=output_file_name))
    else:
        logger.error("delta-filter exited with non-zero code, running for \"{delta_file}\"".format(delta_file=delta_file_name, show_coords_error_log=delta_filter_stderr))
        logger.error("\tdelta-filter error log is stored in \"{delta_filter_error_log}\"".format(delta_file=delta_file_name, delta_filter_error_log=delta_filter_stderr))
    return exitcode


def exit_program():
    logger.critical("An error was encountered and `--ensure-all` flag was set to true. fasta2camsa_points.py exists.")
    exit(1)


def get_assembly_points_from_aligned_contigs(coords_entries, strategy):
    if strategy == "mid-point-sort":
        entries = sorted(coords_entries, key=lambda it: it.mid_coordinate)
        return zip(entries[:-1], entries[1:])
    if strategy == "sliding-window":
        to_go = deque(sorted(coords_entries, key=lambda it: (it.scaffold_start, it.scaffold_end)))
        processing = set()
        finished = []
        pairs = []
        for current in to_go:
            to_be_finished = sorted([entry for entry in processing if entry.scaffold_end <= current.scaffold_start], key=lambda i: i.scaffold_end)
            if len(to_be_finished) == 0 and len(finished) != 0:
                pairs.append((finished[-1], current))
            for item in to_be_finished:
                pairs.append((item, current))
                processing.remove(item)
                finished.append(item)
            processing.add(current)
        return pairs


def filter_fully_covered_contigs(contigs):
    start_entries = [(contig.scaffold_start, contig, "start") if contig.fragment_orientation == "+" else (contig.scaffold_end, contig, "start") for contig in contigs]
    end_entries = [(contig.scaffold_end, contig, "end") if contig.fragment_orientation == "-" else (contig.scaffold_start, contig, "end") for contig in contigs]
    entries = sorted(start_entries + end_entries, key=lambda it: it[0])
    result_chain = []
    processing = deque()
    if "3600" in {contig.fragment_name for contig in contigs}:
        a = 5
    for entry in entries:
        if entry[2] == "start":
            processing.append(entry[1])
        else:
            if entry[1] == processing[0]:
                result_chain.append(entry[1])
                processing.popleft()
            else:
                processing.remove(entry[1])
    return result_chain


if __name__ == "__main__":
    full_description = camsa.full_description_template.format(
        names=camsa.CAMSA_AUTHORS,
        affiliations=camsa.AFFILIATIONS,
        dummy=" ",
        tool="Converting FASTA formatted scaffolding results for further CAMSA processing.",
        information="For more information refer to {docs}".format(docs=camsa.CAMSA_DOCS_URL),
        contact=camsa.CONTACT)
    full_description = "=" * 80 + "\n" + full_description + "=" * 80 + "\n"
    parser = configargparse.ArgParser(description=full_description, formatter_class=configargparse.RawTextHelpFormatter,
                                      default_config_files=[os.path.join(camsa.root_dir, "utils", "fasta", "fasta2camsa_points.ini"),
                                                            os.path.join(camsa.root_dir, "logging.ini")])

    parser.add_argument("contigs", metavar="CONTIGS", help="fasta formatted file with contigs, that served as input for scaffolding purposes")
    parser.add_argument("scaffolds", metavar="SCAFFOLDS", nargs="+", help="fasta formatted result files of contigs scaffolding")
    parser.add_argument("--version", action="version", version=camsa.VERSION)
    parser.add_argument("-c", "--config", is_config_file=True, help="Config file overwriting some of the default settings as well as any flag starting with \"--\".")
    parser.add_argument("-o", "--output-dir", metavar="OUTPUT_DIR", dest="output_dir", default="output", help="Output directory to store temporary and final files.\nDEFAULT: ./output/")
    parser.add_argument("--tmp-dir", metavar="TMP_DIR", dest="tmp_dir", default=None, help="Directory with all intermediate (.delta | .coords) files.")
    parser.add_argument("--overwrite", dest="overwrite", action="store_true", help="Disregards already present \"*.delta\" as well as \"*.coords\" and runs all the stages from scratch.\nDEFAULT: False")
    parser.add_argument("--ensure-all", dest="ensure_all", default=False, type=bool, help="Flag indicating, that is any subprocess of this converter fails, than the whole operation will be canceled.\nDEFAULT: False")

    parser.add_argument("--nucmer-cli-arguments", dest="nucmer_cli_arguments",
                        help="Command line argument, used when running \"nucmer\" software\nDEFAULT: -maxmatch -c 100")

    parser.add_argument("--show-coords-cli-arguments", dest="show_coords_cli_arguments",
                        help="Command line argument, used when running show-coords software\nDEFAULT: -r -c -l")

    parser.add_argument("--delta-filter-cli-arguments", dest="delta_filter_cli_arguments",
                        help="Command line argument, used when running delta-filter software\nDEFAULT: -r -q")

    #######################################################################################################################
    # will definitely make two options this work, once the python bug is fixed: http://bugs.python.org/issue15112
    # by now all the delta files, that match the contigs files name will results in alignment stage for those file to be
    #######################################################################################################################
    # parser.add_argument("--delta-files", default=None, nargs="*", help="paths to the nucmer \"*.delta\" files. For each file present, corresponding alignment stage will be skipped")
    # parser.add_argument("--coords-files", default=None, nargs="*", help="paths to the nucmer \"*.coords\" files. For each file present, corresponding show-coords stage will be skipped")

    parser.add_argument("--nucmer-path", dest="nucmer", help="full path to the nucmer executable.\nDEFAULT: nucmer")
    parser.add_argument("--show-coords-path", dest="show_coords", help="full path to the show-coords executable.\nDEFAULT: show-coords")
    parser.add_argument("--delta-filter-path", dest="delta_filter", help="full path to the delta-filter executable.\nDEFAULT: delta-filter")

    parser.add_argument("--c-cov-threshold", default=90.0, type=float, dest="c_cov_threshold",
                        help="lower coverage bound with respect to each aligned contig. All contigs with coverage less than the threshold are omitted.\nDEAFULT: 90.0")
    parser.add_argument("--c-no-collapse-cons-alignments", default=True, action="store_false", dest="collapse_consecutive_alignments")
    parser.add_argument("--c-keep-fully-covered-contigs", action="store_true", default=False,
                        help="Whether to keep contigs, that are mapped some other contigs, or not.\nDEFAULT: False")
    parser.add_argument("--c-coords-pairs-strategy", choices=["mid-point-sort", "sliding-window", "all"], default="mid-point-sort", dest="coords_to_pairs_strategy",
                        help="A strategy that determines on how assembly pairs from contigs on each scaffold are inferred.\n"
                             "\"mid-point-sort\" -- all contig mapping on each scaffold are sorted by their mid coordinate (start + end) / 2.\n\tSorted sequence of contigs determines n-1 assembly points.\n"
                             "\"sliding-window\" -- all pairs of adjacent extremities of non overlapping contigs will be reported as assembly points."
                             "\nDEFAULT: mid-point-sort ")
    parser.add_argument("--c-logging-level", type=int,
                        choices=[logging.NOTSET, logging.DEBUG, logging.INFO, logging.WARNING, logging.ERROR, logging.CRITICAL],
                        help="Logging level for the converter.\nDEFAULT: {info}".format(info=logging.INFO))
    parser.add_argument("--c-logging-formatter-entry",
                        help="Format string for python logger.")
    args = parser.parse_args()
    start_time = datetime.datetime.now()

    logger = logging.getLogger("CAMSA.utils.fasta2camsa_pairs")
    ch = logging.StreamHandler()
    ch.setLevel(args.c_logging_level)
    logger.setLevel(args.c_logging_level)
    logger.addHandler(ch)
    logger.info(full_description)
    logger.info(parser.format_values())
    ch.setFormatter(logging.Formatter(args.c_logging_formatter_entry))
    logger.info("Starting the converting process")

    args.output_dir = os.path.expanduser(args.output_dir)
    args.output_dir = os.path.abspath(args.output_dir)
    args.tmp_dir = os.path.join(args.output_dir, "fasta2camsa") if args.tmp_dir is None else os.path.abspath(os.path.expanduser(args.tmp_dir))
    args.logs_dir = os.path.join(args.tmp_dir, "logs")

    if not os.path.exists(args.output_dir):
        logger.debug("Output directory \"{directory}\" doesn't exists. Creating one.".format(directory=args.output_dir))
        os.makedirs(args.output_dir)
    if not os.path.exists(args.tmp_dir):
        logger.debug("Creating \"tmp\" directory, where intermediate tools results will be stored")
        os.makedirs(args.tmp_dir)
    if not os.path.exists(args.logs_dir):
        logger.debug("Creating log directory, that would contain all the external tools logs information")
        os.makedirs(args.logs_dir)

    files_in_tmp_dir = [f for f in os.listdir(args.tmp_dir) if isfile(os.path.join(args.tmp_dir, f))]
    delta_files_base_names = {os.path.splitext(f)[0] for f in files_in_tmp_dir if os.path.splitext(f)[1] == ".delta"}
    coords_files_base_names = {os.path.splitext(f)[0] for f in files_in_tmp_dir if os.path.splitext(f)[1] == ".coords"}

    for scaffolds_file in args.scaffolds:
        prefix = get_file_prefix(file_name=scaffolds_file)
        logger.info("Working with \"{scaffolds_file}\"".format(scaffolds_file=scaffolds_file))
        ##################################################################################################
        # obtaining NUCmer alignments results in a form of *.delta file.
        ##################################################################################################
        if (prefix in delta_files_base_names or prefix in coords_files_base_names) and not args.overwrite:
            logger.info("File \"{prefix}.delta\" or \"{prefix}.coords\" is present in the tmp folder. "
                        "Skipping alignment [show-coords] stage for \"{scaffold_file}\" scaffolds file."
                        "".format(prefix=prefix, scaffold_file=scaffolds_file))
        else:
            logger.info("Running NUCmer for \"{scaffolds_file}\" scaffolds file, using \"{contigs_file}\" as query. This might take time."
                        "".format(contigs_file=args.contigs, scaffolds_file=scaffolds_file))
            logger.debug("Results will be stored in \"{prefix}.delta\"."
                         "".format(prefix=prefix))
            result = run_nucmer(contigs_file_name=args.contigs, reference_file_name=scaffolds_file,
                                nucmer_executable_path=args.nucmer,
                                output_dir=args.tmp_dir,
                                logs_dir=args.logs_dir,
                                cli_arguments=args.nucmer_cli_arguments)
            if result != 0 and args.ensure_all:
                exit_program()

        ##################################################################################################
        # obtaining NUCmer alignment results in a form of *.coord file.
        ##################################################################################################
        if prefix in coords_files_base_names and not args.overwrite:
            logger.info("File \"{prefix}.coords\" is present in the tmp folder."
                        " Skipping show-coords stage for \"{prefix}.delta\" file."
                        "".format(prefix=prefix))
        else:
            delta_file_name = [f for f in os.listdir(args.tmp_dir)]
            delta_file_name = [f for f in delta_file_name if os.path.isfile(os.path.join(args.tmp_dir, f))]
            delta_file_name = [f for f in delta_file_name if os.path.splitext(f)[1] == ".delta"]
            delta_file_name = [f for f in delta_file_name if os.path.splitext(f)[0] == prefix]
            if len(delta_file_name) == 0:
                logger.error("Delta file for prefix=\"{prefix}\" was not found in the output folder.".format(prefix=prefix))
                if args.ensure_all:
                    exit_program()
                continue
            original_name_prefix = ".".join(delta_file_name[0].split(".")[:-1])
            delta_file_name = os.path.join(args.tmp_dir, original_name_prefix + ".delta")
            logger.info("Running delta-filter util for \"{delta_file}\" file.".format(delta_file=delta_file_name))
            logger.debug("Results will be stored in \"{prefix}.filtered.delta\"".format(prefix=prefix))
            result = run_delta_filter(delta_file_name=delta_file_name, output_dir=args.tmp_dir, logs_dir=args.logs_dir,
                                      delta_filter_executable_path=args.delta_filter,
                                      cli_arguments=args.delta_filter_cli_arguments)
            if result != 0 and args.ensure_all:
                exit_program()

            delta_file_name = os.path.join(args.tmp_dir, original_name_prefix + ".filtered.delta")
            logger.info("Running show-coords util for \"{delta_file}\" file.".format(delta_file=delta_file_name))
            logger.debug("Results will be stored in \"{prefix}.coords\"".format(prefix=original_name_prefix))
            result = run_show_coords(delta_file_name=delta_file_name, output_dir=args.tmp_dir, logs_dir=args.logs_dir, show_coords_executable_path=args.show_coords, cli_arguments=args.show_coords_cli_arguments)
            if result != 0 and args.ensure_all:
                exit_program()

        ##################################################################################################
        # translating *.coords language of NUCmer into the language of assembly pairs of CAMSA
        # also generates a fragments' length mapping file
        ##################################################################################################
        coords_file = os.path.join(args.tmp_dir, prefix + ".coords")
        if not os.path.exists(coords_file):
            logger.error("Coords file for \"{prefix}\" doesn't exist in the tmp output directory.".format(prefix=prefix))
            if args.ensure_all:
                exit_program()
            continue
        with open(coords_file, "rt") as source:
            logger.info("Parsing \"{coords_file}\".".format(coords_file=coords_file))
            parser = CoordsParser(source=source, lower_frag_cover=args.c_cov_threshold, collapse_consecutive=args.collapse_consecutive_alignments)
            chains = parser.parse_data()

            result_file = os.path.join(args.output_dir, prefix + ".camsa.points")
            with open(result_file, "wt") as destination:
                logger.info("Writing coords data in terms of CAMSA assembly points in \"{camsa_input_file}\"".format(camsa_input_file=result_file))
                writer = csv.writer(destination, delimiter="\t")
                writer.writerow(['origin', 'seq1', 'seq1_or', 'seq2', 'seq2_or', 'gap_size', 'cw'])
                for scaffold, contigs_chain in chains.items():
                    if not args.c_keep_fully_covered_contigs:
                        contigs_chain = filter_fully_covered_contigs(contigs=contigs_chain)
                    assembly_points = get_assembly_points_from_aligned_contigs(coords_entries=contigs_chain,
                                                                               strategy=args.coords_to_pairs_strategy)
                    for left, right in assembly_points:
                        writer.writerow([prefix,
                                         left.fragment_name, left.fragment_orientation,
                                         right.fragment_name, right.fragment_orientation,
                                         min(right.scaffold_start, right.scaffold_end) - max(left.scaffold_end, left.scaffold_start),
                                         "?"])
            logger.info("Finished converting data for \"{prefix}\"".format(prefix=prefix))
    logger.info("Elapsed time: {el_time}".format(el_time=str(datetime.datetime.now() - start_time)))
