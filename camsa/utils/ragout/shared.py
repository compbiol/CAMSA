# -*- coding: utf-8 -*-


def filter_duplications(blocks_by_ids):
    for block_id in list(blocks_by_ids.keys()):
        blocks = blocks_by_ids[block_id]
        genomes = [block.parent_seq.genome_name for block in blocks]
        genomes_set = set(genomes)
        if len(genomes_set) != len(genomes):
            del blocks_by_ids[block_id]


def get_all_genomes_from_blocks(blocks_as_ids):
    genomes = set()
    for block_id, blocks in blocks_as_ids.items():
        for block in blocks:
            genomes.add(block.parent_seq.genome_name)
    return genomes


def filter_indels(blocks_by_ids, all_genomes_as_set=None):
    if all_genomes_as_set is None:
        all_genomes_as_set = get_all_genomes_from_blocks(blocks_as_ids=blocks_by_ids)
    for block_id in list(blocks_by_ids.keys()):
        blocks = blocks_by_ids[block_id]
        genomes_set = {block.parent_seq.genome_name for block in blocks}
        if len(all_genomes_as_set) != len(genomes_set):
            del blocks_by_ids[block_id]


def filter_blocks_by_good_genomes(blocks_by_ids, good_genomes):
    for block_id in list(blocks_by_ids.keys()):
        blocks = blocks_by_ids[block_id]
        new_blocks = [block for block in blocks if block.parent_seq.genome_name in good_genomes]
        if len(new_blocks) == 0:
            del blocks_by_ids[block_id]
        else:
            blocks_by_ids[block_id] = new_blocks


def filter_blocks_by_bad_genomes(blocks_by_ids, bad_genomes):
    for block_id in list(blocks_by_ids.keys()):
        blocks = blocks_by_ids[block_id]
        new_blocks = [block for block in blocks if block.parent_seq.genome_name not in bad_genomes]
        if len(new_blocks) == 0:
            del blocks_by_ids[block_id]
        else:
            blocks_by_ids[block_id] = new_blocks
