#!/usr/bin/python
# -*- coding: UTF-8 -*-

import argparse

import pandas as pd

from orthomap import datasets, gtf2t2g

file = datasets.mouse_gtf(datapath="/tmp")
expected_columns = [
    "gene_id",
    "gene_id_version",
    "transcript_id",
    "transcript_id_version",
    "gene_name",
    "gene_type",
    "protein_id",
    "protein_id_version",
]


def test_define_parse():
    parse = gtf2t2g.define_parser()
    assert isinstance(parse, argparse.ArgumentParser)


def test_gtf2t2g_only_filename():
    output = gtf2t2g.parse_gtf(file, q=True)
    assert (output.columns == expected_columns).all()
    assert isinstance(output, pd.DataFrame)
    empty_columns = [
        "gene_id_version",
        "transcript_id_version",
        "gene_name",
        "gene_type",
        "protein_id",
        "protein_id_version",
    ]
    for col in empty_columns:
        assert len(output[col].unique()) == 1
        assert output[col].unique()[0] == "None" or output[col].unique()[0] is None


def test_gtf2t2g_with_gene_type_and_versions():
    output = gtf2t2g.parse_gtf(file, g=True, b=True, v=True, q=True)
    assert (output.columns == expected_columns).all()
    assert isinstance(output, pd.DataFrame)
    empty_columns = ["protein_id", "protein_id_version"]
    for col in empty_columns:
        assert len(output[col].unique()) == 1
        assert output[col].unique()[0] is None


def test_gtf2t2g_with_gene_protein_id_version():
    output = gtf2t2g.parse_gtf(file, g=True, p=True, q=True)
    assert (output.columns == expected_columns).all()
    assert isinstance(output, pd.DataFrame)
    empty_columns = [
        "gene_id_version",
        "transcript_id_version",
        "gene_type",
        "protein_id_version",
    ]
    for col in empty_columns:
        assert len(output[col].unique()) == 1
        assert output[col].unique()[0] == "None" or output[col].unique()[0] is None


def test_mus_musculus():
    expected_first_row = [
        "ENSMUSG00000000001",
        "ENSMUSG00000000001.5",
        "ENSMUST00000000001",
        "ENSMUST00000000001.5",
        "Gnai3",
        None,
        "ENSMUSP00000000001",
        "ENSMUSP00000000001.5",
    ]
    expected_last_row = [
        "ENSMUSG00002076992",
        "ENSMUSG00002076992.1",
        "ENSMUST00020181762",
        "ENSMUST00020181762.1",
        "7SK",
        None,
        None,
        None,
    ]
    output = gtf2t2g.parse_gtf(file, g=True, p=True, s=True, q=True, v=True)
    assert (output.iloc[0].values == expected_first_row).all()
    assert (output.iloc[-1].values == expected_last_row).all()
