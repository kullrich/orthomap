.. _get_orhtomap:

Get query species taxonomic lineage information
===============================================

Given a species name or taxonomic ID, the query species lineage information is
extracted with the help of the `ete3` python toolkit and the `NCBI taxonomy`
([Huerta-Cepas et al., 2016](https://doi.org/10.1093/molbev/msw046)).

This information is needed alongside with the taxonomic classifications for all
species used in the OrthoFinder comparison.

The `orthomap` submodule `qlin` helps to get this information for you with the
`qlin.get_qlin()` function as follows::

    >>> from orthomap import qlin
    >>> qlin.get_qlin(q = 'Danio rerio')
    >>> qlin.get_qlin(qt = '7955')


Get query species orthomap from OrthoFinder results
===================================================