#! /usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import print_function

import argparse
import datetime
import itertools
import logging
import os
import shutil
import sys
from collections import defaultdict

import six
from jinja2 import Template

sys.path.append(os.path.dirname(__file__))
sys.path.append(os.path.dirname(os.path.dirname(__file__)))

import camsa
from core import io as camsa_io
from core import merging
from core.comparative_analysis import compute_and_update_assembly_points_conflicts
from core.data_structures import Assembly, assign_ids_to_assembly_points, merge_assembly_points, assign_parents_to_children
from core.merging import MergingStrategies, update_assembly_points_with_merged_assembly


if __name__ == "__main__":
    full_description = camsa.full_description_template.format(
        names=camsa.CAMSA_AUTHORS,
        affiliations=camsa.AFFILIATIONS,
        dummy=" ",
        tool="CAMSA is a tool for Comparative Analysis and Merging of Scaffold Assemblies",
        information="For more information refer to wiki at github.com/aganezov/camsa/wiki",
        contact=camsa.CONTACT)
    full_description = "=" * 80 + "\n" + full_description + "=" * 80 + "\n"

    parser = argparse.ArgumentParser(description=full_description, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("input", nargs="+")
    parser.add_argument("--c-default-exact-cw", dest="cw_exact", type=float, default=1.0)
    parser.add_argument("--c-default-candidate-cw", dest="cw_candidate", type=float, default=0.9)
    parser.add_argument("--c-min-cw-threshold", dest="min_cw", type=float, default=0.0)
    parser.add_argument("--c-merging-strategy", choices=[MergingStrategies.progressive_merging.value, MergingStrategies.maximal_matching.value],
                        default=MergingStrategies.maximal_matching.value, dest="merging_strategy")
    parser.add_argument("--c-merging-allow-cycles", action="store_false", dest="acyclic", default=True)
    parser.add_argument("--version", action="version", version=camsa.VERSION)
    parser.add_argument("-o", "--c-output-dir", dest="output_report_dir", default=None)
    parser.add_argument("--c-logging-level", dest="logging_level", default=logging.INFO, type=int,
                        choices=[logging.NOTSET, logging.DEBUG, logging.INFO, logging.WARNING, logging.ERROR, logging.CRITICAL],
                        help="Logging level for the converter.\nDEFAULT: {info}".format(info=logging.INFO))
    args = parser.parse_args()

    start_time = datetime.datetime.now()

    #######################################
    #           logging setup             #
    #######################################
    logger = logging.getLogger("CAMSA.main")
    ch = logging.StreamHandler()
    logger.setLevel(args.logging_level)
    logger.addHandler(ch)
    logger.info(full_description)
    ch.setFormatter(camsa.formatter)
    logger.info("Starting the analysis")

    #######################################
    #           input stage               #
    #######################################
    logger.info("Processing input")

    # key is the origin of assembly point; values is the list of all points from that source
    assembly_points_by_sources = camsa_io.read_assembly_points_from_input_sources(sources=args.input,
                                                                                  default_cw_eae=args.cw_exact,
                                                                                  default_cw_cae=args.cw_candidate)
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
    #        comparative analysis         #
    #######################################
    logger.info("Computing assembly points conflicts")
    compute_and_update_assembly_points_conflicts(assembly_points_by_ids=merged_assembly_points_by_ids)

    #######################################
    #       merging assemblies            #
    #######################################
    logger.info("Obtaining a merged assembly, using {strategy} strategy".format(strategy=args.merging_strategy))
    merged_assembly_graph = merging.strategies_bindings[args.merging_strategy](assembly_points_by_sources=assembly_points_by_sources,
                                                                               acyclic=args.acyclic,
                                                                               min_cw=args.min_cw)
    update_assembly_points_with_merged_assembly(original_assembly_points_by_ids=original_assembly_points_by_ids,
                                                merged_assembly_points_by_ids=merged_assembly_points_by_ids,
                                                merged_assembly_graph=merged_assembly_graph)

    #######################################
    #           output stage              #
    #######################################
    logger.info("Preparing output")
    if args.output_report_dir is None:
        args.output_report_dir = os.path.join(os.getcwd(), "camsa_{date}".format(
            date=datetime.datetime.now().strftime("%b_%d_%Y__%H_%M")))
        logger.debug("Output directory was not specified. Automatically generated name: \"{out_dir}\"."
                     "".format(out_dir=args.output_report_dir))

    args.output_report_dir = os.path.abspath(os.path.expanduser(args.output_report_dir))
    if not os.path.exists(args.output_report_dir):
        os.makedirs(args.output_report_dir)

    if os.path.exists(os.path.join(args.output_report_dir, "libs")) and os.path.isdir(
            os.path.join(args.output_report_dir, "libs")):
        shutil.rmtree(os.path.join(args.output_report_dir, "libs"))
    root_dir = os.path.dirname(os.path.realpath(__file__))
    shutil.copytree(src=os.path.join(root_dir, "libs"), dst=os.path.join(args.output_report_dir, "libs"))
    output_html_report_file_name = os.path.join(args.output_report_dir, "report.html")

    template = Template(source=open(os.path.join(root_dir, "report_template.html"), "rt").read())

    individual_assemblies.sort(key=lambda it: it.name.lower())
    assemblies_to_ids = {assembly.name: cnt for cnt, assembly in enumerate(individual_assemblies, start=1)}
    sources = [assembly.name for assembly in individual_assemblies]

    assemblies_to_colors = {assemblies_to_ids[source]: color for source, color in zip(sources, ['red', 'blue', 'green', 'purple',
                                                                                                'orange', 'pink', 'brown', 'navy', 'steelblue'])}

    with open(output_html_report_file_name, "wt") as dest:
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
                    "lengths": []
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
                }
            }), file=dest)
    logger.info("CAMSA report is written to \"{output_report_file}\"".format(output_report_file=output_html_report_file_name))
    logger.info("Finished Comparative Analysis and Merging of input assemblies.")
    logger.info("Elapsed time: {el_time}".format(el_time=str(datetime.datetime.now() - start_time)))
    logger.info("Thank you for using CAMSA!")
