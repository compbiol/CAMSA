# -*- coding: utf-8 -*-
import os
import re

import sys
from collections import defaultdict

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))))

from camsa.utils.ragout.data_structures import RagoutSequence, Block

block_id_pattern = re.compile("Block #(?P<block_id>\d+)")
ENTRIES_SEPARATOR = "--------------------------------------------------------------------------------"


def read_from_file(path, silent_fail=False, delimiter="\t"):
    with open(path, "rt") as source:
        whole_text = source.read()
        sequences_by_ids, blocks_by_ids = parse_ragout_coords(ragout_coords_as_a_string=whole_text, delimiter=delimiter, silent_fail=silent_fail)
        return sequences_by_ids, blocks_by_ids


def parse_ragout_coords(ragout_coords_as_a_string, delimiter="\t", silent_fail=False):
    data = [string_entry.strip() for string_entry in ragout_coords_as_a_string.split(ENTRIES_SEPARATOR) if len(string_entry.strip()) > 0]
    seq_part = data[0]
    blocks = data[1:]
    sequences_by_ids = parse_ragout_seq_part(entry_as_a_string=seq_part, delimiter=delimiter, silent_fail=silent_fail)
    all_blocks = []
    for block_string in blocks:
        current_blocks = parse_ragout_block_entry(entry_as_a_string=block_string, sequences_by_ids=sequences_by_ids, silent_fail=silent_fail, delimiter=delimiter)
        all_blocks.extend(current_blocks)
    blocks_by_ids = defaultdict(list)
    for block in all_blocks:
        blocks_by_ids[block.name].append(block)
    return sequences_by_ids, blocks_by_ids


def get_seq_headers_columns_indexes(headers_string, delimiter, silent_fail=False):
    """
            Seq_id	Size	Description
    """
    headers_string = headers_string.strip()
    data = headers_string.split(delimiter)
    seq_id_columns = data.index("Seq_id")
    length_column = data.index("Size")
    description_column = data.index("Description")
    if any(map(lambda x: x < 0, [seq_id_columns, length_column, description_column])) and not silent_fail:
        raise ValueError("Could not parse block headers")
    return {"seq_id_column": seq_id_columns, "seq_length_column": length_column, "seq_description_column": description_column}


def parse_ragout_seq_part(entry_as_a_string, delimiter="\t", silent_fail=False):
    entries = [entry.strip() for entry in entry_as_a_string.split("\n")]
    seq_headers = entries[0]
    seq_entries = entries[1:]
    columns_indexes = get_seq_headers_columns_indexes(headers_string=seq_headers, delimiter=delimiter, silent_fail=silent_fail)
    if any(map(lambda x: x < 0, columns_indexes.values())) and silent_fail:
        return
    result_sequences = {}
    for entry_string in seq_entries:
        seq = parse_ragout_seq_entry(entry_as_a_string=entry_string,
                                     seq_id_column=columns_indexes["seq_id_column"],
                                     seq_length_column=columns_indexes["seq_length_column"],
                                     seq_description_column=columns_indexes["seq_description_column"], delimiter=delimiter)
        result_sequences[seq.ragout_id] = seq
    return result_sequences


def parse_ragout_seq_entry(entry_as_a_string, seq_id_column=0, seq_length_column=1, seq_description_column=2, delimiter="\t"):
    """
        Seq_id	Size	Description
        27	11234	Anc0.Anc0refChr7043
    """
    data = entry_as_a_string.split(delimiter)
    seq_id = int(float(data[seq_id_column]))
    seq_length = int(float(data[seq_length_column]))
    seq_name = data[seq_description_column]
    return RagoutSequence(name=seq_name, ragout_id=seq_id, length=seq_length)


def parse_ragout_block_entry(entry_as_a_string, sequences_by_ids=None, silent_fail=False, delimiter="\t"):
    """
        Block #4873
        Seq_id	Strand	Start	End	Length
        18386	+	36450191	36456417	6226
        19238	-	6907	699	6208
        24172	-	9270	2915	6355
        28652	-	606707	600419	6288
        30888	-	606707	600419	6288
        ....
    """
    entries = [entry.strip() for entry in entry_as_a_string.split("\n")]
    block_id_string = entries[0]
    columns_headers = entries[1]
    blocks = entries[2:]
    block_id = parse_block_id_string(data_string=block_id_string, silent_fail=silent_fail)
    if block_id is None and silent_fail:
        return
    columns_indexes = get_blocks_headers_columns_indexes(headers_string=columns_headers, delimiter=delimiter, silent_fail=silent_fail)
    if any(map(lambda x: x < 0, columns_indexes.values())) and silent_fail:
        return
    result_blocks = []
    for block_string in blocks:
        result_blocks.append(parse_ragout_block_string(data_string=block_string, block_id=block_id,
                                                       seq_id_column=columns_indexes["seq_id_column"],
                                                       strand_column=columns_indexes["strand_column"],
                                                       start_column=columns_indexes["start_column"],
                                                       end_column=columns_indexes["end_column"]))
    if sequences_by_ids is not None:
        for block in result_blocks:
            try:
                block.parent_seq = sequences_by_ids[block.parent_seq]
            except KeyError:
                if not silent_fail:
                    raise
    return result_blocks


def get_blocks_headers_columns_indexes(headers_string, delimiter="\t", silent_fail=False):
    """
        Seq_id	Strand	Start	End	Length
    """
    headers_string = headers_string.strip()
    data = headers_string.split(delimiter)
    seq_id_columns = data.index("Seq_id")
    strand_column = data.index("Strand")
    start_column = data.index("Start")
    end_column = data.index("End")
    if any(map(lambda x: x < 0, [seq_id_columns, strand_column, start_column, end_column])) and not silent_fail:
        raise ValueError("Could not parse block headers")
    return {"seq_id_column": seq_id_columns, "strand_column": strand_column, "start_column": start_column, "end_column": end_column}


def parse_block_id_string(data_string, pattern=block_id_pattern, silent_fail=False):
    """
        Block #xxxx
    """
    data_string = data_string.strip()
    match = pattern.match(data_string)
    if match is None:
        if not silent_fail:
            raise ValueError("Block id string was not successfully matched by the regex pattern")
    else:
        block_id = int(match.group("block_id"))
        return block_id


def parse_ragout_block_string(data_string, block_id, seq_id_column=0, strand_column=1, start_column=2, end_column=3, delimiter="\t"):
    """
        Seq_id	Strand	Start	End	Length
        18386	+	36450191	36456417	6226
    """
    data_string = data_string.strip()
    data = data_string.split(delimiter)
    seq_id = int(data[seq_id_column])
    strand = data[strand_column]
    start = int(data[start_column])
    end = int(data[end_column])
    if strand == "-":
        start, end = end, start
    block = Block(name=block_id, start=start, end=end, strand=strand, parent_seq=seq_id)
    return block
