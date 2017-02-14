#! /usr/bin/env python
# -*- coding: utf-8 -*-
import datetime
import enum
import logging
import os
import sys
from collections import defaultdict

import configargparse

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))))

import camsa
import camsa.core.io as camsa_io
from camsa.core.data_structures import AssemblyPoint


class Component(object):
    def __init__(self, object_id, object_beg, object_end, part_number, component_type):
        self.object_id = object_id
        self.object_beg = object_beg
        self.object_end = object_end
        self.part_number = part_number
        self.component_type = component_type

    @classmethod
    def from_agp_data(cls, data):
        raise NotImplementedError("Only specific heirs have implementation for this method")

    def is_gap_component(self):
        raise NotImplementedError("Only specific heirs have implementation for this method")

    @staticmethod
    def is_scaffold_component(component_type):
        return component_type not in ["U", "N"]


class ScaffoldComponent(Component):
    def __init__(self, object_id, object_beg, object_end, part_number, component_type, component_id, component_beg, component_end, orientation):
        super(ScaffoldComponent, self).__init__(object_id=object_id, object_beg=object_beg, object_end=object_end, part_number=part_number, component_type=component_type)
        self.component_id = component_id
        self.component_beg = component_beg
        self.component_end = component_end
        self.orientation = orientation

    @classmethod
    def from_agp_data(cls, data):
        object_id = data[0]
        object_beg = int(data[1])
        object_end = int(data[2])
        part_number = int(data[3])
        component_type = data[4]
        component_id = data[5]
        component_beg = int(data[6])
        component_end = int(data[7])
        orientation = data[8]
        return cls(object_id=object_id, object_beg=object_beg, object_end=object_end, part_number=part_number, component_type=component_type,
                   component_id=component_id, component_beg=component_beg, component_end=component_end, orientation=orientation)

    @property
    def is_gap_component(self):
        return False


class GapComponent(Component):
    def __init__(self, object_id, object_beg, object_end, part_number, component_type, gap_length, gap_type, linkage, linkage_evidence):
        super(GapComponent, self).__init__(object_id=object_id, object_beg=object_beg, object_end=object_end, part_number=part_number, component_type=component_type)
        self.gap_length = gap_length
        self.gap_type = gap_type
        self.linkage = linkage
        self.linkage_evidence = linkage_evidence

    @classmethod
    def from_agp_data(cls, data):
        object_id = data[0]
        object_beg = int(data[1])
        object_end = int(data[2])
        part_number = int(data[3])
        component_type = data[4]
        gap_length = int(data[5])
        gap_type = data[6]
        linkage = data[7]
        linkage_evidence = data[8]
        return cls(object_id=object_id, object_beg=object_beg, object_end=object_end, part_number=part_number, component_type=component_type,
                   gap_length=gap_length, gap_type=gap_type, linkage=linkage, linkage_evidence=linkage_evidence)

    @property
    def is_gap_component(self):
        return True


def get_index_of_first_non_gap_element(component_sequence):
    for cnt, component in enumerate(component_sequence):
        if not component.is_gap_component:
            return cnt
    return -1


def get_next_non_gap_element(components, current_index):
    for i in range(current_index + 1, len(components)):
        next = components[i]
        if not next.is_gap_component:
            return i, next
    return None, None


def get_gap_size(components, current_index, next_index):
    if current_index == next_index - 1:
        return 1
    gap_components = components[current_index+1:next_index]
    result = 0
    for component in gap_components:
        assert component.is_gap_component
        if component.component_type == "U":
            return "?"
        result += component.gap_length
    return result


def camsa_orientation_from_agp(agp_or):
    if agp_or in ["?", "0", "na"]:
        return "?"
    return agp_or


def object_as_camsa_points(object_as_components, extra_data):
    result = []
    components = object_as_components
    f_non_gap_element = get_index_of_first_non_gap_element(component_sequence=components)
    current_element = components[f_non_gap_element]
    current_index = f_non_gap_element
    next_index, next_element = get_next_non_gap_element(components=components, current_index=current_index)
    while next_element is not None:
        gap_size = get_gap_size(components, current_index, next_index)
        seq1_or = camsa_orientation_from_agp(agp_or=current_element.orientation)
        seq2_or = camsa_orientation_from_agp(agp_or=next_element.orientation)
        assembly_point = AssemblyPoint(seq1=current_element.component_id, seq2=next_element.component_id,
                                       seq1_or=seq1_or, seq2_or=seq2_or, sources=[extra_data.origin],
                                       gap_size=gap_size, cw="?")
        result.append(assembly_point)
        current_element = next_element
        current_index = next_index
        next_index, next_element = get_next_non_gap_element(components=components, current_index=current_index)
    return result


if __name__ == "__main__":
    full_description = camsa.full_description_template.format(
        names=camsa.CAMSA_AUTHORS,
        affiliations=camsa.AFFILIATIONS,
        dummy=" ",
        tool="Converting AGPv2.0 formatted scaffolding results for further CAMSA processing.",
        information="For more information refer to {docs}".format(docs=camsa.CAMSA_DOCS_URL),
        contact=camsa.CONTACT)
    full_description = "=" * 80 + "\n" + full_description + "=" * 80 + "\n"
    parser = configargparse.ArgParser(description=full_description, formatter_class=configargparse.RawTextHelpFormatter,
                                      default_config_files=[os.path.join(camsa.root_dir, "utils", "agp", "agp2camsa_points.ini"),
                                                            os.path.join(camsa.root_dir, "logging.ini")])
    parser.add_argument("-c", "--config", is_config_file=True,
                        help="Config file overwriting some of the default settings as well as any flag starting with \"--\".")
    parser.add_argument("--version", action="version", version=camsa.VERSION)
    parser.add_argument("--c-logging-level", type=int,
                        choices=[logging.NOTSET, logging.DEBUG, logging.INFO, logging.WARNING, logging.ERROR, logging.CRITICAL],
                        help="Logging level for the converter.\nDEFAULT: {info}".format(info=logging.INFO))
    parser.add_argument("--c-logging-formatter-entry",
                        help="Format string for python logger.")

    parser.add_argument("agp", type=configargparse.FileType("rt"), nargs="+", default=sys.stdin,
                        help="Input stream of AGPv2 formatted scaffold assemblies\nDEFAULT: stdin")
    parser.add_argument("--origin", type=str, default=None,
                        help="Identifier for the assembly, that would be specified in the \"origin\" column of CAMSA points\nDEFAULT: inferred from input file names")
    parser.add_argument("--o-format", type=str,
                        help="The CAMSA-out formatting for the assembly points obtained form the AGPv2 formatted scaffold assemblies")
    parser.add_argument("-o", "--output", type=configargparse.FileType("wt"), default=sys.stdout,
                        help="The stream where CAMSA formatted assembly points are outputted\nDEFAULT: stdout")

    args = parser.parse_args()

    start_time = datetime.datetime.now()
    #######################################
    #           logging setup             #
    #######################################
    logger = logging.getLogger("CAMSA.utils.agp2camsa_points")
    ch = logging.StreamHandler()
    ch.setLevel(args.c_logging_level)
    logger.setLevel(args.c_logging_level)
    logger.addHandler(ch)
    logger.info(full_description)
    logger.info(parser.format_values())
    ch.setFormatter(logging.Formatter(args.c_logging_formatter_entry))
    logger.info("Starting the converting process")

    objects = defaultdict(list)
    assembly_points = []

    if args.origin is None:
        logger.debug("\"origin\" were not specified explicitly for this AGP data. Inferring from the data stream")
        sources = []
        for source in args.agp:
            sources.append(os.path.splitext(os.path.basename(str(source.name)))[0])
        args.origin = ".".join(sources)
        logger.debug("\"origin\" has been inferred to \"{sources}\"".format(sources=args.origin))

    #######################################
    #      reading input AGP data         #
    #######################################
    logger.info("Reading input AGP formatted data")
    for source in args.agp:
        for line in source:
            line = line.strip()
            if line.startswith("#"):
                logger.debug("Skipping comment line: {line}".format(line=line))
                continue
            data = line.split("\t", 8)
            if Component.is_scaffold_component(data[4]):
                logger.debug("Processing a non-gap data line: {line}".format(line=line))
                component = ScaffoldComponent.from_agp_data(data=data)
            else:
                logger.debug("Processing a gap data line: {line}".format(line=line))
                component = GapComponent.from_agp_data(data=data)
            objects[component.object_id].append(component)

    for object_id in objects.keys():
        logger.debug("Processing object {object_id}".format(object_id=object_id))
        objects[object_id] = sorted(objects[object_id], key=lambda component: (component.object_beg, component.object_end))
        components = objects[object_id]
        camsa_points = object_as_camsa_points(object_as_components=components, extra_data=args)
        assembly_points.extend(camsa_points)

    logger.info("Writing CAMSA formatted assembly poitns to {file}".format(file=args.output.name))
    camsa_io.write_assembly_points(assembly_points=assembly_points, destination=args.output, output_setup=args.o_format)
    logger.info("Finished the conversion.")
    logger.info("Elapsed time: {el_time}".format(el_time=str(datetime.datetime.now() - start_time)))
