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
    parser.add_argument("--c-logging-level", type=int,
                        choices=[logging.NOTSET, logging.DEBUG, logging.INFO, logging.WARNING, logging.ERROR, logging.CRITICAL],
                        help="Logging level for the converter.\nDEFAULT: {info}".format(info=logging.INFO))
    parser.add_argument("--c-logging-formatter-entry",
                        help="Format string for python logger.")

    parser.add_argument("agp", type=configargparse.FileType("rt"), default=sys.stdin)
    parser.add_argument("-o", "--output", type=configargparse.FileType("wt"), default=sys.stdout)

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

    #######################################
    #      reading input AGP data         #
    #######################################
    logger.info("Reading input AGP formatted data")
    for line in args.agp:
        line = line.strip()
        if line.startswith("#"):
            continue    # skipping comment lines
        data = line.split("\t", 8)
        if Component.is_scaffold_component(data[4]):
            component = ScaffoldComponent.from_agp_data(data=data)
        else:
            component = GapComponent.from_agp_data(data=data)
        objects[component.object_id].append(component)

    logger.info("Elapsed time: {el_time}".format(el_time=str(datetime.datetime.now() - start_time)))
