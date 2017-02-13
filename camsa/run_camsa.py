#! /usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import print_function

import datetime
import itertools
import logging
import os
import shutil
import sys
from collections import defaultdict

import configargparse
import six
from jinja2 import FileSystemLoader
from jinja2.environment import Environment

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import camsa
from camsa.core import io as camsa_io
from camsa.core import merging
from camsa.core.reference_analysis import analyze_and_update_assembly_points_based_on_reference
from camsa.core.comparative_analysis import compute_and_update_assembly_points_conflicts
from camsa.core.data_structures import Assembly, assign_ids_to_assembly_points, merge_assembly_points, assign_parents_to_children, to_json, Sequence
from camsa.core.merging import MergingStrategies, update_assembly_points_with_merged_assembly, update_gap_sizes_in_merged_assembly

if __name__ == "__main__":
    full_description = camsa.full_description_template.format(
        names=camsa.CAMSA_AUTHORS,
        affiliations=camsa.AFFILIATIONS,
        dummy=" ",
        tool="CAMSA is a tool for Comparative Analysis and Merging of Scaffold Assemblies",
        information="For more information refer to {docs}".format(docs=camsa.CAMSA_DOCS_URL),
        contact=camsa.CONTACT)
    full_description = "=" * 80 + "\n" + full_description + "=" * 80 + "\n"

    parser = configargparse.ArgParser(description=full_description,
                                      formatter_class=configargparse.RawTextHelpFormatter,
                                      default_config_files=[os.path.join(camsa.root_dir, "run_camsa.ini"),
                                                            os.path.join(camsa.root_dir, "logging.ini")])
    parser.add_argument("points", nargs="+",
                        help="A list of input files, representing a standard CAMSA format for assembly points.")
    parser.add_argument("--seqi", default=None,
                        help="A file with sequences' information about fragments, involved in input assembly points")
    parser.add_argument("--seqi-delimiter", type=str, default="\t",
                        help="A delimiter character for the file containing sequences' information.\nDEFAULT: \\t")
    parser.add_argument("--seqi-ensure-all", action="store_true", dest="lengths_ensure", default=False,
                        help="Ensures, that either all, or non of the sequences participating in assembly points have information provided for them.\nDEFAULT: False")
    parser.add_argument("-c", "--config", is_config_file=True,
                        help="Config file path with settings for CAMSA to run with.\nOverwrites the default CAMSA configuration file.\nValues in config file can be overwritten by command line arguments.")
    parser.add_argument("--c-cw-exact", type=float,
                        help="A confidence weight value assigned to oriented assembly points and respective exact assembly edges,\nin case \"?\" is specified as the respective assembly point confidence weight.\nDEFAULT: 1.0")
    parser.add_argument("--c-cw-candidate", type=float,
                        help="A confidence weight value assigned to semi/un-oriented assembly points and respective candidate assembly edges,\nin case \"?\" is specified as the respective assembly point confidence weight.\nDEFAULT: 0.75")
    # parser.add_argument("--c-comparative-disable", action="store_false", dest="comparative", default=True,
    #                     help="")
    parser.add_argument("--ref-disable", action="store_false", dest="reference", default=True,
                        help="")
    parser.add_argument("--ref", type=str, default="", dest="reference_name",
                        help="")
    # parser.add_argument("--c-merging-disable", action="store_false", dest="merging", default=True,
    #                     help="")
    parser.add_argument("--c-merging-cw-min", type=float,
                        help="A threshold for the minimum cumulative confidence weight for merged assembly edges in MSAG.\nEdges with confidence weight below are not considered in the \"merged\" assembly construction.\nDEFAULT: 0.0")
    parser.add_argument("--c-merging-strategy", choices=[MergingStrategies.greedy_merging.value, MergingStrategies.maximal_matching.value],
                        default=MergingStrategies.maximal_matching.value,
                        help="A strategy to produced a merged assembly from the given ones.\nDEFAULT: maximal-matching")
    parser.add_argument("--c-merging-cycles", dest="allow_cycles", action="store_true", default=False,
                        help="Whether to allow cycles in the produced merged assembly.\nDEFAULT: False")
    parser.add_argument("--version", action="version", version=camsa.VERSION)
    parser.add_argument("--i-delimiter", default="\t", type=str,
                        help="String used as a delimiter in the input files with CAMSA assembly points")
    parser.add_argument("-o", "--o-dir",
                        help="A directory, where CAMSA will store all of the produced output (report, assets, etc).\nDEFAULT: camsa_{date}")
    # parser.add_argument("--o-interactive-disable", action="store_false", dest="o_interactive", default=True,
    #                     help="")
    parser.add_argument("--o-merged-format", type=str,
                        help="The CAMSA-out formatting for the merged scaffold assemblies in a form of CAMSA points.")
    parser.add_argument("--o-subgroups-format", type=str,
                        help="The CAMSA-out formatting for the merged scaffold assemblies in a form of CAMSA points.")
    parser.add_argument("--o-collapsed-format", type=str,
                        help="The CAMSA-out formatting for the collapsed assembly points and their computed conflicts.")
    parser.add_argument("--o-original-format", type=str,
                        help="The CAMSA-out formatting for the non-collapsed assembly points and their computed conflicts.")
    parser.add_argument("--c-logging-level", default=logging.INFO, type=int,
                        choices=[logging.NOTSET, logging.DEBUG, logging.INFO, logging.WARNING, logging.ERROR, logging.CRITICAL],
                        help="Logging level for CAMSA.\nDEFAULT: {info}".format(info=logging.INFO))
    parser.add_argument("--c-logging-formatter-entry",
                        help="Format string for python logger.")
    args = parser.parse_args()

    start_time = datetime.datetime.now()

    #######################################
    #           logging setup             #
    #######################################
    if args.o_dir is None:
        args.o_dir = os.path.join(os.getcwd(), "camsa_{date}".format(
            date=datetime.datetime.now().strftime("%b_%d_%Y__%H_%M")))

    args.output_dir = os.path.abspath(os.path.expanduser(args.o_dir))
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    logger = logging.getLogger("CAMSA.main")
    ch = logging.StreamHandler()
    logger.setLevel(args.c_logging_level)
    logger.addHandler(ch)
    fh = logging.FileHandler(os.path.join(args.output_dir, "camsa.log"), "wt")
    fh.setLevel(logging.DEBUG)
    logger.addHandler(fh)
    logger.info(full_description)
    logger.info(parser.format_values())
    ch.setFormatter(logging.Formatter(args.c_logging_formatter_entry))
    fh.setFormatter(logging.Formatter(args.c_logging_formatter_entry))
    logger.info("Starting the analysis")

    #######################################
    #           input stage               #
    #######################################
    logger.info("Processing input")

    # key is the origin of assembly point; values is the list of all points from that source
    logger.info("Reading assembly points")
    assembly_points_by_sources = camsa_io.read_assembly_points_from_input_sources(sources=args.points,
                                                                                  default_cw_eae=args.c_cw_exact,
                                                                                  default_cw_cae=args.c_cw_candidate,
                                                                                  delimiter=args.i_delimiter)
    or_seqi = defaultdict(list)
    if args.seqi is not None:
        args.seqi = os.path.abspath(os.path.expanduser(args.seqi))
        logger.info("Reading sequences' info")
        with open(args.seqi, "rt") as source:
            camsa_io.read_seqi_from_input_sources(source=source, delimiter=args.seqi_delimiter, destination=or_seqi)
    seqi = {}
    for seq_name in list(or_seqi.keys()):
        entry_list = or_seqi[seq_name]
        seqi[seq_name] = entry_list[0]
    if len(seqi) > 0:
        sequences_ids_from_assembly_points = set()
        for source, aps in assembly_points_by_sources.items():
            for ap in aps:
                sequences_ids_from_assembly_points.add(ap.seq1)
                sequences_ids_from_assembly_points.add(ap.seq2)
        sequences_ids_from_lengths = set(seqi.keys())
        difference = sequences_ids_from_assembly_points - sequences_ids_from_lengths
        if len(difference) > 0:
            logger.warning("Some sequences, that participate in assembly points do not have lengths associated with them,"
                           "while others do.")
            logger.warning("Problematic sequences: {seqs}".format(seqs=",".join(sorted(difference))))
            if args.seqi_ensure:
                logger.error("A flag \"--seqi-ensure-all\" was set, so the program terminates")
            else:
                logger.warning("Assigning lengths of -1 to all problematic sequences")
                for seq_id in difference:
                    seqi[seq_id] = Sequence(name=seq_id, length=-1)

    #######################################
    #       extracting reference          #
    #######################################
    reference_assembly = None
    if args.reference and args.reference_name != "":
        try:
            reference_aps = assembly_points_by_sources.pop(args.reference_name)
            reference_assembly = Assembly(name=args.reference_name, aps=reference_aps)
        except KeyError:
            logger.critical("Supplied reference \"{reference_name}\" was not found among assembly sources [{avail_sources}]".format(reference_name=args.reference_name, avail_sources=",".join(assembly_points_by_sources.keys())))
            exit(1)

    id_generator = itertools.count()
    original_assembly_points_by_ids = {}
    for assembly_points in assembly_points_by_sources.values():
        original_assembly_points_by_ids.update(assign_ids_to_assembly_points(assembly_points=assembly_points, id_prefix="or_",
                                                                             id_generator=id_generator))

    #######################################
    #       assembly points merging       #
    #######################################
    logger.info("Merging assembly points from different sources into a set of unique ones.")
    merged_assembly_points = merge_assembly_points(assembly_points_by_source=assembly_points_by_sources)
    merged_assembly_points_by_ids = assign_ids_to_assembly_points(assembly_points=merged_assembly_points, id_prefix="m_")
    assign_parents_to_children(children_assembly_points_by_ids=original_assembly_points_by_ids,
                               parent_assembly_points_by_ids=merged_assembly_points_by_ids)

    #######################################
    #       assemblies subgroups          #
    #######################################
    logger.info("Processing assemblies' subgroups")
    tmp_individual_assemblies = defaultdict(list)
    for ap in merged_assembly_points:
        for source_name in ap.sources:
            tmp_individual_assemblies[source_name].append(ap)
    individual_assemblies = [Assembly(name=name, aps=aps) for name, aps in tmp_individual_assemblies.items()]

    grouped_assemblies = []
    for i in range(1, len(individual_assemblies) + 1):
        for assembly_combination in itertools.combinations(sorted([a.name for a in individual_assemblies]), i):
            assembly = Assembly(name=assembly_combination, aps=[])
            for ap in merged_assembly_points:
                if len(ap.sources) == len(assembly_combination) and all(map(lambda entry: entry in ap.sources, assembly_combination)):
                    assembly.aps.append(ap)
            grouped_assemblies.append(assembly)
    grouped_assemblies = list(filter(lambda a: len(a.aps) > 0, sorted(grouped_assemblies, key=lambda entry: len(entry.aps), reverse=True)))

    #######################################
    #        reference   analysis         #
    #######################################
    # for assembly in individual_assemblies:
    #     analyze_and_update_assembly_points_based_on_reference(assembly=assembly, ref_assembly=reference_assembly)

    #######################################
    #        comparative analysis         #
    #######################################
    logger.info("Computing assembly points conflicts")
    compute_and_update_assembly_points_conflicts(assembly_points_by_ids=merged_assembly_points_by_ids)

    #######################################
    #       merging assemblies            #
    #######################################
    logger.info("Obtaining a merged assembly, using {strategy} strategy".format(strategy=args.c_merging_strategy))
    merged_assembly_graph = merging.strategies_bindings[args.c_merging_strategy](assembly_points_by_sources=assembly_points_by_sources,
                                                                                 acyclic=not args.allow_cycles,
                                                                                 min_cw=args.c_merging_cw_min)
    update_assembly_points_with_merged_assembly(original_assembly_points_by_ids=original_assembly_points_by_ids,
                                                merged_assembly_points_by_ids=merged_assembly_points_by_ids,
                                                merged_assembly_graph=merged_assembly_graph)
    update_gap_sizes_in_merged_assembly(original_assembly_points_by_ids=original_assembly_points_by_ids,
                                        merged_assembly_points_by_ids=merged_assembly_points_by_ids)

    #######################################
    #           output stage              #
    #######################################
    logger.info("Preparing output")

    # copying assets required for the HTML report
    libs_report_dir = os.path.join(args.output_dir, "libs")
    camsa_io.remove_dir(dir_path=libs_report_dir)
    shutil.copytree(src=os.path.join(camsa.root_dir, "libs"), dst=libs_report_dir)
    output_html_report_file_name = os.path.join(args.output_dir, "report.html")

    # "input" subdir of the report
    # will contain a configuration as well as assembly points files
    input_report_dir = os.path.join(args.output_dir, "input")
    camsa_io.remove_dir(dir_path=input_report_dir)
    os.makedirs(input_report_dir)
    input_report_config_path = os.path.join(input_report_dir, "camsa_config.txt")
    with open(input_report_config_path, "wt") as destination:
        print("# NOTE: this is not a valid config, but rather a summary of the utilized options", file=destination)
        print(parser.format_values(), file=destination)
    for pairs_path in args.points:
        full_path = os.path.abspath(os.path.expanduser(pairs_path))
        base_name = os.path.basename(full_path)
        shutil.copyfile(src=full_path, dst=os.path.join(input_report_dir, base_name))

    # "merged" subdir of the report
    # will contain assembly points, that constitute the merged assembly
    merged_report_dir = os.path.join(args.output_dir, "merged")
    camsa_io.remove_dir(dir_path=merged_report_dir)
    os.makedirs(merged_report_dir)
    merged_report_points_path = os.path.join(merged_report_dir, "merged.camsa.points")
    with open(merged_report_points_path, "wt") as destination:
        camsa_io.write_assembly_points(destination=destination,
                                       assembly_points=[ap for ap in merged_assembly_points if ap.participates_in_merged],
                                       output_setup=args.o_merged_format)

    # "comparative" subdir of the report
    # will contain assembly points divided into subgroups based in the agreement in input assemblies
    comparative_report_dir = os.path.join(args.output_dir, "comparative")
    camsa_io.remove_dir(dir_path=comparative_report_dir)
    subgroups_report_dir = os.path.join(comparative_report_dir, "subgroups")
    camsa_io.remove_dir(comparative_report_dir)
    os.makedirs(comparative_report_dir)
    os.makedirs(subgroups_report_dir)
    for group in grouped_assemblies:
        comparative_report_group_points_path = os.path.join(subgroups_report_dir, "{group_name}.camsa.points".format(group_name=".".join(group.name)))
        with open(comparative_report_group_points_path, "wt") as destination:
            camsa_io.write_assembly_points(assembly_points=group.aps,
                                           destination=destination,
                                           output_setup=args.o_subgroups_format)

    original_points_path = os.path.join(comparative_report_dir, "original.camsa.points")
    with open(original_points_path, "wt") as destination:
        camsa_io.write_assembly_points(assembly_points=original_assembly_points_by_ids.values(),
                                       destination=destination,
                                       output_setup=args.o_original_format)

    collapsed_points_path = os.path.join(comparative_report_dir, "collapsed.camsa.points")
    with open(collapsed_points_path, "wt") as destination:
        camsa_io.write_assembly_points(destination=destination,
                                       assembly_points=merged_assembly_points,
                                       output_setup=args.o_collapsed_format)
    env = Environment()
    env.filters['tojson'] = to_json

    individual_assemblies.sort(key=lambda it: it.name.lower())
    assemblies_to_ids = {assembly.name: "A" + str(cnt) for cnt, assembly in enumerate(individual_assemblies, start=1)}
    sources = [assembly.name for assembly in individual_assemblies]

    assemblies_to_colors = {assemblies_to_ids[source]: color for source, color in zip(sources, ['red', 'blue', 'green', 'purple',
                                                                                                'orange', 'pink', 'brown', 'navy', 'steelblue'])}

    with open(output_html_report_file_name, "wt") as dest:
        env.loader = FileSystemLoader(os.path.join(camsa.root_dir, "html"))
        template = env.get_template("base_template.html")
        six.print_(template.render(
            data={
                "assemblies": individual_assemblies,
                "assemblies_intersections": [],
                "assemblies_conflicts": [],
                "graph_compiled": False,
                "aps": merged_assembly_points,
                "assemblies_to_ids": assemblies_to_ids,
                "assemblies_to_colors": assemblies_to_colors,
                "grouped_assemblies": grouped_assemblies,
                "fragments": {
                    "seqi": seqi,
                    'max_length': max([seq.length for seq in seqi.values()]) if len(seqi) > 0 else -1,
                    'min_length': min([seq.length for seq in seqi.values()]) if len(seqi) > 0 else -1,
                }
            },
            settings={
                "cytoscape": {
                    "draw_timeout": 300000
                }
            },
            meta={
                "camsa": {
                    "version": camsa.VERSION
                },
                "date": start_time
            }), file=dest)
    logger.info("CAMSA report is written to \"{output_report_file}\"".format(output_report_file=output_html_report_file_name))
    logger.info("Finished Comparative Analysis and Merging of input assemblies.")
    end_time = datetime.datetime.now()
    logger.info("Elapsed time: {el_time}".format(el_time=str(end_time - start_time)))
    logger.info("Thank you for using CAMSA!")
