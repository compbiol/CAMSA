# -*- coding: utf-8 -*-
import csv
import os
from collections import defaultdict

from camsa.core.data_structures import AssemblyPoint


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


def read_pairs(source, delimiter="\t", destination=None, default_cw_eae=1, default_cw_cae=0.9):
    """

    :param source: file like object tot APs data from
    :param delimiter: tab/comma/etc separator
    :param destination: data structure, where information about APs will be stored
    :param default_cw_eae: confidence weight for exact AE, in case ? is provided in source
    :param default_cw_cae: confidence wight for candidate AE, in case ? is provided in source
    :return: destination data structure, that can be viewed as a default dict of list of APs, where key is the source of the AP
    """
    def extract_nullable_numerical_value(field, row, fn_relations):
        value = row.get(fn_relations[field], "?")
        try:
            value = float(value)
        except ValueError:
            value = "?"
        return value

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
    fn_relations = get_fn_relations_for_column_names(fieldnames=fieldnames, aliases=PAIRS_COLUMN_ALIASES)
    for row in filter(lambda entry: not entry[fieldnames[0]].startswith("#"), reader):
        seq_id = row(fn_relations["seq_id"])
        seq_length = row(fn_relations["seq_length"])
        destination[seq_id] = seq_length
    return destination
