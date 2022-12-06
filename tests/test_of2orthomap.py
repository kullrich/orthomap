import pandas as pd

from orthomap import of2orthomap

def test_of2orthomap_continuity_False():
    query_orthomap, orthofinder_species_list, of_species_abundance = of2orthomap.get_orthomap(
    seqname='Danio_rerio.GRCz11.cds.longest',
    qt='7955',
    sl='ensembl_105_orthofinder_species_list.tsv',
    oc='ensembl_105_orthofinder_Orthogroups.GeneCount.tsv',
    og='ensembl_105_orthofinder_Orthogroups.tsv',
    continuity=False,
    quite=True)

    assert isinstance(query_orthomap, pd.DataFrame)
    assert (query_orthomap.columns == ['seqID', 'Orthogroup', 'PSnum', 'PStaxID', 'PSname']).all()
    assert isinstance(orthofinder_species_list, pd.DataFrame)
    assert isinstance(of_species_abundance, pd.DataFrame)

    # Ask KU about the relationship of these three dataframes.
    # check the index in the `of_species_abundance`

def test_of2orthomap_continuity_True():
    query_orthomap, orthofinder_species_list, of_species_abundance = of2orthomap.get_orthomap(
    seqname='Danio_rerio.GRCz11.cds.longest',
    qt='7955',
    sl='ensembl_105_orthofinder_species_list.tsv',
    oc='ensembl_105_orthofinder_Orthogroups.GeneCount.tsv',
    og='ensembl_105_orthofinder_Orthogroups.tsv',
    continuity=True,
    quite=True)

    assert isinstance(query_orthomap, pd.DataFrame)
    assert (query_orthomap.columns == ['seqID', 'Orthogroup', 'PSnum', 'PStaxID', 'PSname', 'PScontinuity']).all()
    assert isinstance(orthofinder_species_list, pd.DataFrame)
    assert isinstance(of_species_abundance, pd.DataFrame)

def test_orthomap_last_entrance():
    # Is it true that the last column should always be the queried species?
    pass