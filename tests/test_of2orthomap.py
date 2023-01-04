#!/usr/bin/python
# -*- coding: UTF-8 -*-

import argparse

import pandas as pd

from orthomap import datasets, of2orthomap

ensembl105_oc, ensembl105_og, ensembl105_sl = datasets.ensembl105('/tmp')


def test_define_parse():
    parse = of2orthomap.define_parser()
    assert isinstance(parse, argparse.ArgumentParser)

def test_of2orthomap_continuity_false():
    query_orthomap, orthofinder_species_list, of_species_abundance = of2orthomap.get_orthomap(
        seqname='Danio_rerio.GRCz11.cds.longest',
        qt='7955',
        oc=ensembl105_oc,
        og=ensembl105_og,
        sl=ensembl105_sl,
        continuity=False,
        quite=True)
    assert isinstance(query_orthomap, pd.DataFrame)
    assert (query_orthomap.columns == ['seqID', 'Orthogroup', 'PSnum', 'PStaxID', 'PSname']).all()
    assert isinstance(orthofinder_species_list, pd.DataFrame)
    assert isinstance(of_species_abundance, pd.DataFrame)
    # Ask KU about the relationship of these three dataframes.
    # check the index in the `of_species_abundance`


def test_of2orthomap_continuity_true():
    query_orthomap, orthofinder_species_list, of_species_abundance = of2orthomap.get_orthomap(
        seqname='Danio_rerio.GRCz11.cds.longest',
        qt='7955',
        oc=ensembl105_oc,
        og=ensembl105_og,
        sl=ensembl105_sl,
        continuity=True,
        quite=True)
    assert isinstance(query_orthomap, pd.DataFrame)
    assert (query_orthomap.columns == ['seqID', 'Orthogroup', 'PSnum', 'PStaxID', 'PSname', 'PScontinuity']).all()
    assert isinstance(orthofinder_species_list, pd.DataFrame)
    assert isinstance(of_species_abundance, pd.DataFrame)
