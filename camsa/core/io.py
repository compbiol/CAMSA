# -*- coding: utf-8 -*-
import csv
import os
import shutil
from collections import defaultdict
from enum import Enum

from camsa.core.data_structures import AssemblyPoint, APFieldOutExtractorConverter, Sequence


def get_fn_relations_for_column_names(fieldnames, aliases):
    canonical_field_names = [aliases.get(name.lower(), name) for name in fieldnames]
    fn_relations = {value: key for key, value in aliases.items()}
    fn_relations.update({canonical: name for name, canonical in zip(fieldnames, canonical_field_names)})
    return fn_relations


PAIRS_COLUMN_ALIASES = {
    ########################
    "species": "origin",
    "organisms": "origin",
    ########################
    "ctg1": "seq1",
    ########################
    "ctg2": "seq2",
    ########################
    "orientation_ctg1": "seq1_or",
    "ctg1_orientation": "seq1_or",
    "ctg1_or": "seq1_or",
    ########################
    "ctg2_orientation": "seq2_or",
    "orientation_ctg2": "seq2_or",
    "ctg2_or": "seq2_or",
    ########################
    "distance": "gap_size",
    "ctg1-ctg2_gap": "gap_size",
    ########################
    "score": "cw",
}


def extract_nullable_value(field, row, fn_relations, default="?"):
    if field not in fn_relations:
        return default
    return row.get(fn_relations[field], default)


def extract_nullable_numerical_value(field, row, fn_relations, default="?"):
    value = extract_nullable_value(field=field, row=row, fn_relations=fn_relations, default=default)
    try:
        value = float(value)
    except (ValueError, TypeError):
        value = default
    return value


def read_pairs(source, delimiter="\t", destination=None, default_cw_eae=1, default_cw_cae=0.9):
    """

    :param source: file like object tot APs data from
    :param delimiter: tab/comma/etc separator
    :param destination: data structure, where information about APs will be stored
    :param default_cw_eae: confidence weight for exact AE, in case ? is provided in source
    :param default_cw_cae: confidence wight for candidate AE, in case ? is provided in source
    :return: destination data structure, that can be viewed as a default dict of list of APs, where key is the source of the AP
    """


    if destination is None:
        destination = defaultdict(list)
    reader = csv.DictReader(source, delimiter=delimiter)
    fieldnames = reader.fieldnames
    fn_relations = get_fn_relations_for_column_names(fieldnames=fieldnames, aliases=PAIRS_COLUMN_ALIASES)
    for row in filter(lambda entry: not entry[fieldnames[0]].startswith("#"), reader):
        origin = row[fn_relations["origin"]]
        cw = extract_nullable_numerical_value(field="cw", row=row, fn_relations=fn_relations)
        if cw == "?":
            cw = default_cw_eae if "?" not in [row[fn_relations["seq1_or"]], row[fn_relations["seq2_or"]]] else default_cw_cae
        gap_size = extract_nullable_numerical_value(field="gap_size", row=row, fn_relations=fn_relations)
        seq1 = row[fn_relations["seq1"]]
        seq2 = row[fn_relations["seq2"]]
        if seq1 == seq2:
            # no support for duplicated seqs present
            continue
        destination[origin].append(AssemblyPoint(seq1=seq1,
                                                 seq2=seq2,
                                                 seq1_or=row[fn_relations["seq1_or"]],
                                                 seq2_or=row[fn_relations["seq2_or"]],
                                                 sources=[row[fn_relations["origin"]]],
                                                 cw=cw,
                                                 gap_size=gap_size))
    return destination


def read_assembly_points_from_input_sources(sources, delimiter="\t", default_cw_eae=1, default_cw_cae=0.75):
    """

    :param sources: list of file paths with input AP data
    :param delimiter: tab/comma/etc separator
    :param default_cw_eae: confidence weight for exact AE, in case ? is provided in source
    :param default_cw_cae: confidence wight for candidate AE, in case ? is provided in source
    :return: destination data structure, that can be viewed as a default dict of list of APs, where key is the source of the AP
    """
    result = defaultdict(list)
    for file_name in sources:
        file_name = os.path.abspath(os.path.expanduser(file_name))
        with open(file_name, "rt") as source:
            read_pairs(source=source, delimiter=delimiter, destination=result,
                       default_cw_eae=default_cw_eae, default_cw_cae=default_cw_cae)
    return result


LENGTHS_COLUMN_ALIASES = {
    ########################
    "ctg_id": "seq_id",
    ########################
    "ctg_length": "seq_length"
}


def read_lengths(source, delimiter="\t", destination=None):
    if destination is None:
        destination = {}
    reader = csv.DictReader(source, delimiter=delimiter)
    fieldnames = reader.fieldnames
    fn_relations = get_fn_relations_for_column_names(fieldnames=fieldnames, aliases=LENGTHS_COLUMN_ALIASES)
    for row in filter(lambda entry: not entry[fieldnames[0]].startswith("#"), reader):
        seq_id = row[fn_relations["seq_id"]]
        seq_length = int(float(row[fn_relations["seq_length"]]))
        destination[seq_id] = Sequence(name=seq_id, length=seq_length)
    return destination


class OrientationChoice(Enum):
    original = 0
    merged = 1


def get_header_and_extract_list(settings):
    data = settings.split("|")
    converters_setups = []
    for entry in data:
        converters_setups.append(entry.split(","))
    header = [entry[0] for entry in converters_setups]
    converters = [APFieldOutExtractorConverter(field_name=entry[1], converter_name=entry[2]) for entry in converters_setups]
    return header, converters


def write_assembly_points(assembly_points, destination, output_setup, delimiter="\t"):
    """ Output a collection of assembly point in a text format to the specified stream

    :param assembly_points: an iterable with a collection of assembly points to be written down
    :param destination: a file like object to write to
    :param delimiter: a separator used for the SV format
    :param orientation_type: a choice for AP relative seq orientations to be displayed (original = input vs inferred = merged).
        Makes a difference only for the un/semi-oriented APs
    """
    intra_delimiter = "," if delimiter != "," else ";"

    writer = csv.writer(destination, delimiter=delimiter)
    header, field_processor = get_header_and_extract_list(settings=output_setup)
    writer.writerow(header)
    for ap in assembly_points:
        assembly_points_entries_list = [processor.extract_field_value_str(ap) for processor in field_processor]
        writer.writerow(assembly_points_entries_list)


def write_seqi(sequences, destination, output_setup, delimiter="\t"):
    writer = csv.writer(destination, delimiter=delimiter)
    header, field_processor = get_header_and_extract_list(settings=output_setup)
    writer.writerow(header)
    for seq in sequences:
        seq_entries_list = [processor.extract_field_value_str(seq) for processor in field_processor]
        writer.writerow(seq_entries_list)


def read_seqi_from_input_sources(source, delimiter="\t", destination=None):
    if destination is None:
        destination = defaultdict(list)
    reader = csv.DictReader(source, delimiter=delimiter)
    fieldnames = reader.fieldnames
    fn_relations = get_fn_relations_for_column_names(fieldnames=fieldnames, aliases=LENGTHS_COLUMN_ALIASES)
    for row in filter(lambda entry: not entry[fieldnames[0]].startswith("#"), reader):
        seq_id = row[fn_relations["seq_id"]]
        length = extract_nullable_numerical_value(field="length", row=row, fn_relations=fn_relations, default=None)
        parent_seq_id = extract_nullable_value(field="parent_seq_id", row=row, fn_relations=fn_relations, default=None)
        start = extract_nullable_numerical_value(field="start", row=row, fn_relations=fn_relations, default=None)
        end = extract_nullable_numerical_value(field="end", row=row, fn_relations=fn_relations, default=None)
        strand = extract_nullable_value(field="strand", row=row, fn_relations=fn_relations, default="+")
        annotation = extract_nullable_value(field="annotation", row=row, fn_relations=fn_relations, default=None)
        seq_group_id = extract_nullable_value(field="seq_group_id", row=row, fn_relations=fn_relations, default="Default")
        seq = Sequence(name=seq_id, length=length, parent_seq_id=parent_seq_id, start=start, end=end, strand=strand, annotation=annotation, seq_group_id=seq_group_id)
        destination[seq.name].append(seq)
    return destination


def remove_dir(dir_path):
    if os.path.exists(dir_path) and os.path.isdir(dir_path):
        shutil.rmtree(dir_path)
