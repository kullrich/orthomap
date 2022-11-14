import numpy as np
import pandas as pd

from orthomap import gtf2t2g

file = "Mus_musculus.GRCm39.108.chr.gtf.gz"

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
        assert output[col].unique()[0] == "None" or output[col].unique()[0] == None


def test_gtf2t2g_with_gene_type_and_versions():
    output = gtf2t2g.parse_gtf(file, g=True, b=True, v=True, q=True)

    assert (output.columns == expected_columns).all()
    assert isinstance(output, pd.DataFrame)

    empty_columns = ["protein_id", "protein_id_version"]

    for col in empty_columns:
        assert len(output[col].unique()) == 1
        assert output[col].unique()[0] == None


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
        assert output[col].unique()[0] == "None" or output[col].unique()[0] == None
