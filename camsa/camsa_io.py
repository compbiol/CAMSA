# -*- coding: utf-8 -*-
import csv
from collections import defaultdict

from camsa.data_structures import AssemblyPoint

aliases = {
    ########################
    "species": "origin",
    "organisms": "origin",
    ########################
    "ctg1": "ctg1",
    ########################
    "ctg2": "ctg2",
    ########################
    "orientation_ctg1": "ctg1_or",
    "ctg1_orientation": "ctg1_or",
    ########################
    "ctg2_orientation": "ctg2_or",
    "orientation_ctg2": "ctg2_or",
    ########################
    "distance": "gap_size",
    "ctg1-ctg2_gap": "gap_size",
    ########################
    "score": "cw",

}


def read_pairs(source, delimiter="\t", destination=None, default_cw_eae=1, default_cw_pae=0.9):
    if destination is None:
        destination = defaultdict(list)
    reader = csv.DictReader(source, delimiter=delimiter)
    fieldnames = reader.fieldnames
    canonical_field_names = [aliases.get(name.lower(), name) for name in fieldnames]
    fn_relations = {value: key for key, value in aliases.items()}
    fn_relations.update({canonical: name for name, canonical in zip(fieldnames, canonical_field_names)})
    for row in filter(lambda entry: not entry[fieldnames[0]].startswith("#"), reader):
        origin = row[fn_relations["origin"]]
        cw = row.get(fn_relations["cw"], "?")
        try:
            cw = float(cw)
        except ValueError:
            cw = "?"
        if cw == "?":
            cw = default_cw_eae if "?" not in [row[fn_relations["ctg1_or"]], row[fn_relations["ctg2_or"]]] else default_cw_pae
        destination[origin].append(AssemblyPoint(ctg1=row[fn_relations["ctg1"]],
                                                 ctg2=row[fn_relations["ctg2"]],
                                                 ctg1_or=row[fn_relations["ctg1_or"]],
                                                 ctg2_or=row[fn_relations["ctg2_or"]],
                                                 sources=[row[fn_relations["origin"]]],
                                                 cw=cw))
    return destination
