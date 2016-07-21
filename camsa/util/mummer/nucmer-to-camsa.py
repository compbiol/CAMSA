#! /usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse
import logging
import subprocess
import os
from os.path import isfile


def get_file_prefix(file_name):
    return os.path.basename(os.path.splitext(file_name)[0])


def run_nucmer(contigs_file_name, reference_file_name, nucmer_executable_path, cli_arguments, output_dir, nucmer_stdout=None, nucmer_stderr=None):
    prefix = get_file_prefix(reference_file_name)
    if nucmer_stdout is None:
        nucmer_stdout = os.path.join(output_dir, "logs", "nucmer_{prefix}.stdout.txt".format(prefix=prefix))
    if nucmer_stderr is None:
        nucmer_stderr = os.path.join(output_dir, "logs", "nucmer_{prefix}.stderr.txt".format(prefix=prefix))
    cli_args = [nucmer_executable_path] + cli_arguments.split(" ") + ['-p', os.path.join(output_dir, prefix), reference_file_name, contigs_file_name, '>', nucmer_stdout, '2>', nucmer_stderr]
    logger.info("Running:\n\t" + " ".join(cli_args))
    exitcode = subprocess.call(" ".join(cli_args), shell=True)
    if exitcode == 0:
        logger.info("NUCmer finished running for \"{scaffolds_file}\" scaffolds file.\n"
                    "NUCmer alignment output is stored in \"{delta_output}\" file.\n"
                    "NUCmer log is stored in \"{nucmer_log}\".\n---\n"
                    "".format(scaffolds_file=reference_file_name, delta_output=os.path.join(output_dir, "{prefix}.delta".format(prefix=prefix)), nucmer_log=nucmer_stdout))
    else:
        logger.error("NUCmer exited with non-zero code, running for \"{scaffolds}\" scaffolds file. "
                     "NUCmer logs are stored in:\n\tstdout: \"{nucmer_log}\"\n\tstderr: \"{nucmer_err}\"\n---\n".format(scaffolds=reference_file_name, nucmer_log=nucmer_stdout, nucmer_err=nucmer_stderr))
    return exitcode


def run_show_coords(delta_file_name, output_dir, show_coords_executable_path, cli_arguments, show_coords_stderr=None):
    prefix = get_file_prefix(delta_file_name)
    output_file_name = os.path.join(output_dir, prefix + ".coords")
    if show_coords_stderr is None:
        show_coords_stderr = os.path.join(output_dir, "logs", "show_coords_{prefix}.stderr.txt".format(prefix=prefix))
    cli_args = [show_coords_executable_path] + cli_arguments.split(" ") + [delta_file_name, '>', output_file_name, '2>', show_coords_stderr]
    logger.info("Running:\n\t" + " ".join(cli_args))
    exitcode = subprocess.call(" ".join(cli_args), shell=True)
    if exitcode == 0:
        logger.info("show-coords util has finished running for \"{delta_file}\".\n"
                    "show-coords output is stored in \"{coords_output}\" file.\n---\n"
                    "".format(delta_file=delta_file_name, coords_output=output_file_name))
    else:
        logger.error("show-coords exited with non-zero code, running for \"{delta_file}\".\n"
                     "show-coords error log is stored in \"{show_coords_error_log}\"\n---\n".format(delta_file=delta_file_name, show_coords_error_log=show_coords_stderr))
    return exitcode


def exit_program():
    logger.critical("An error was encountered and `--ensure-all` flag was set to true. NUCmer-to-CAMSA exists.")
    exit(1)


if __name__ == "__main__":
    full_description = "-" * 80 + \
                       "\nSergey Aganezov & Max A. Alekseyev (c)\n" + \
                       "Computational Biology Institute, The George Washington University.\n\n" + \
                       "Preparation of fasta formatted scaffolding results for further CAMSA processing.\n" + \
                       "With any questions, please, contact Sergey Aganezov [aganezov(at)gwu.edu].\n" + \
                       "-" * 80
    parser = argparse.ArgumentParser(description=full_description, formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument("contigs", metavar="CONTIGS", help="fasta formatted file with contigs, that served as input for scaffolding purposes")
    parser.add_argument("scaffolds", metavar="SCAFFOLDS", nargs="+", help="fasta formatted result files of contigs scaffolding")
    parser.add_argument("-o", metavar="OUTPUT_DIR", dest="output_dir", default="output", help="output directory to store temporary and final files.\nDEFAULT: ./output/")
    parser.add_argument("--overwrite", dest="overwrite", default=False, type=bool, help="Disregards already present \"*.delta\" as well as \"*.coords\" and runs all the stages from scratch.\nDEFAULT: False")
    parser.add_argument("--ensure-all", dest="ensure_all", default=False, type=bool, help="Flag indicating, that is any subprocess of this converter fails, than the whole operation will be canceled.\nDEFAULT: False")

    parser.add_argument("--nucmer-cli-arguments", dest="nucmer_cli_arguments", default="-maxmatch -c 100")

    parser.add_argument("--show-coords-arguments", dest="show_coords_cli_arguments", default="-r -c -l")

    #######################################################################################################################
    # will definitely make two options this work, once the python bug is fixed: http://bugs.python.org/issue15112
    # by now all the delta files, that match the contigs files name will results in alignment stage for those file to be
    #######################################################################################################################
    # parser.add_argument("--delta-files", default=None, nargs="*", help="paths to the nucmer \"*.delta\" files. For each file present, corresponding alignment stage will be skipped")
    # parser.add_argument("--coords-files", default=None, nargs="*", help="paths to the nucmer \"*.coords\" files. For each file present, corresponding show-coords stage will be skipped")

    parser.add_argument("--nucmer-path", default="/usr/local/bin/nucmer", dest="nucmer", help="full path to the nucmer executable.\nDEFAULT: nucmer")
    parser.add_argument("--show-coords-path", default="/usr/local/bin/show-coords", dest="show_coords", help="full path to the show-coords executable.\nDEFAULT: show-coords")

    parser.add_argument("--c-cov-threshold", default=90.0, type=float, help="lower coverage bound with respect to each aligned contig. All contigs with coverage less than the threshold are omitted.\nDEAFULT: 90.0")

    parser.add_argument("--logging-level", dest="logging_level", default=logging.INFO, choices=[logging.NOTSET, logging.DEBUG, logging.INFO, logging.WARNING, logging.ERROR, logging.CRITICAL],
                        help="Logging level for the converter.\nDEFAULT: {info}".format(info=logging.INFO))
    args = parser.parse_args()

    logger = logging.getLogger("NUCmer-to-CAMSA")
    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)
    logger.setLevel(logging.DEBUG)
    logger.addHandler(ch)
    logger.info("Starting the converting process")

    args.output_dir = os.path.expanduser(args.output_dir)
    args.tmp_dir = os.path.join(args.output_dir, "tmp")
    args.logs_dir = os.path.join(args.output_dir, "logs")

    if not os.path.exists(args.output_dir):
        logger.debug("Output directory \"{directory}\" doesn't exists. Creating one.".format(directory=args.output_dir))
        os.makedirs(args.output_dir)
    if not os.path.exists(os.path.join(args.output_dir), "tmp"):
        logger.debug("Creating \"tmp\" directory, where intermediate tools results will be stored")
        os.makedirs(args.tmp_dir)
    if not os.path.exists(os.path.join(args.output_dir, "logs")):
        logger.debug("Creating log directory, that would contain all the external tools logs information")
        os.makedirs(args.logs_dir)

    files_in_output_dir = [f for f in os.listdir(args.output_dir) if isfile(os.path.join(args.output_dir, f))]
    delta_files_base_names = {os.path.splitext(f)[0] for f in files_in_output_dir if os.path.splitext(f)[1] == ".delta"}
    coords_files_base_names = {os.path.splitext(f)[0] for f in files_in_output_dir if os.path.splitext(f)[1] == ".coords"}

    for scaffolds_file in args.scaffolds:
        prefix = get_file_prefix(file_name=scaffolds_file)
        ##################################################################################################
        # obtaining NUCmer alignments results in a form of *.delta file.
        ##################################################################################################
        if prefix in delta_files_base_names and not args.overwrite:
            logger.info("---\nFile \"{prefix}.delta\" is present in the output folder. Skipping alignment stage for \"{scaffold_file}\" contigs file.\n---\n"
                        "".format(prefix=prefix, scaffold_file=scaffolds_file))
        else:
            logger.info("---\nRunning NUCmer for \"{contigs_file}\" contigs file, using \"{scaffold_file}\" as reference."
                        " Results will be stored in \"{prefix}.delta\".\nThis might take time.\n"
                        "".format(contigs_file=args.contigs, prefix=prefix, scaffold_file=scaffolds_file))
            result = run_nucmer(contigs_file_name=args.contigs, reference_file_name=scaffolds_file,
                                nucmer_executable_path=args.nucmer,
                                output_dir=args.output_dir,
                                cli_arguments=args.nucmer_cli_arguments)
            if result != 0 and args.ensure_all:
                exit_program()

        ##################################################################################################
        # obtaining NUCmer alignment results in a form of *.coord file.
        ##################################################################################################
        if prefix in coords_files_base_names and not args.overwrite:
            logger.info("---\nFile \"{prefix}.coords\" is present in the output folder. Skipping show-coords stage for \"{prefix}\".delta file.\n---\n"
                        "".format(prefix=prefix))
        else:
            delta_file_name = [f for f in os.listdir(args.output_dir)]  # if isfile(os.path.join(args.output_dir, f)) and os.path.splitext(f)[1] == ".delta" and os.path.splitext(f)[0] == prefix]
            delta_file_name = [f for f in delta_file_name if os.path.isfile(os.path.join(args.output_dir, f))]
            delta_file_name = [f for f in delta_file_name if os.path.splitext(f)[1] == ".delta"]
            delta_file_name = [f for f in delta_file_name if os.path.splitext(f)[0] == prefix]
            if len(delta_file_name) == 0:
                logger.error("---\nDelta file for prefix=\"{prefix}\" was not found in the output folder.\n---\n".format(prefix=prefix))
                if args.ensure_all:
                    exit_program()
                continue
            delta_file_name = os.path.join(args.output_dir, delta_file_name[0])
            logger.info("---\nRunning show-coords util for \"{delta_file}\" file."
                        " Results will be stored in \"{prefix}.coords\"".format(prefix=prefix, delta_file=delta_file_name))
            result = run_show_coords(delta_file_name=delta_file_name, output_dir=args.output_dir, show_coords_executable_path=args.show_coords, cli_arguments=args.show_coords_cli_arguments)
            if result != 0 and args.ensure_all:
                exit_program()
