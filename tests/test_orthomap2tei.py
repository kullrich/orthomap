#!/usr/bin/python
# -*- coding: UTF-8 -*-

import pandas as pd

from orthomap import datasets, orthomap2tei


def test_read_orthomap():
    expected_columns = ['GeneID', 'Phylostratum']
    query_orthomap = orthomap2tei.read_orthomap('examples/Sun2021_Orthomap.tsv')
    assert isinstance(query_orthomap, pd.DataFrame)
    assert all(query_orthomap.columns == expected_columns)


def test_geneset_overlap():
    g1 = ['A', 'B', 'C']
    g2 = ['A', 'D', 'E', 'F']
    overlap = orthomap2tei.geneset_overlap(g1, g2)
    assert overlap['g1_g2_overlap'].values == 1
    assert overlap['g1_ratio'].values == 1 / 3
    assert overlap['g2_ratio'].values == 1 / 4


def test_replace_by():
    # Note that with this function you always have to replace all the elements.
    x = ['A', 'B', 'C']
    excepted_x = ['D', 'K', 'E']
    xmatch = ['A', 'C', 'B']
    xreplace = ['D', 'E', 'K']
    newx = orthomap2tei.replace_by(x, xmatch, xreplace)
    assert newx == excepted_x


def test_keep_min_max_case_max():
    # Not sure we need a function for this. It can be done in two lines.
    my_orthomap = pd.DataFrame.from_dict(
        {'GeneID': ['g1', 'g1', 'g2', 'g3', 'g3'], 'Phylostrata': [3, 1, 2, 5, 7]}
    )
    ascending = {'max': False, 'min': True}
    expected_orthomap = (
        my_orthomap.sort_values('Phylostrata', ascending=ascending['max'])
        .drop_duplicates('GeneID')
        .sort_values('Phylostrata', ascending=ascending['max'])
    )
    filtered_orthomap = orthomap2tei.keep_min_max(my_orthomap, keep='max')
    assert (expected_orthomap.values == filtered_orthomap.values).all()


def test_keep_min_max_case_min():
    # Not sure we need a function for this. It can be done in two lines.
    my_orthomap = pd.DataFrame.from_dict(
        {'GeneID': ['g1', 'g1', 'g2', 'g3', 'g3'], 'Phylostrata': [3, 1, 2, 5, 7]}
    )
    ascending = {'max': False, 'min': True}
    expected_orthomap = (
        my_orthomap.sort_values('Phylostrata', ascending=ascending['min'])
        .drop_duplicates('GeneID')
        .sort_values('Phylostrata', ascending=ascending['min'])
    )
    filtered_orthomap = orthomap2tei.keep_min_max(my_orthomap, keep='min')
    assert (expected_orthomap.values == filtered_orthomap.values).all()


def test_split_gene_id_by_gene_age():
    # Not sure what the function is doing.
    pass


def test_get_psd():
    # No data to run the example.
    pass

def test_min_max_to_01():
    ndarray_example_one = [1 for _ in range(5)]
    assert orthomap2tei.min_max_to_01(ndarray_example_one) == [0 for _ in range(5)]
    ndarray_example_two = [1 for _ in range(5)] + [2]
    assert orthomap2tei.min_max_to_01(ndarray_example_two) == [0 for _ in range(5)] + [1]


def test_get_rematrix():
    pass
