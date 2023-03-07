Welcome to orthomap's documentation!
====================================
|GitHub Workflow Status| |PyPI| |PyPI - Python Version| |PyPI - Wheel| |Licence| |ReadTheDocs|

`orthomap` is a python package that can:

- extract orthologous maps from `OrthoFinder <https://github.com/davidemms/OrthoFinder>`_ results.

- extract orthologous maps from `eggnog_6.0 <http://eggnog6.embl.de/download/eggnog_6.0/>`_ results.

- show species lineage information from the `NCBI taxonomy <https://www.ncbi.nlm.nih.gov/taxonomy>`_ database.

- extract gene and transcript information from GTF files (`GTF File Format <https://www.ensembl.org/info/website/upload/gff.html>`_).

- match gene and transcript IDs with scRNA data sets.

- calculate and visualize transcriptome evolutionary index (TEI) on scRNA data sets.

- calculate and visualize partial transcriptome evolutionary index (TEI) profiles on scRNA data sets.

- calculate mean/relative expression per gene age class.

- show number of expressed genes per gene age class.

Source code is available at `orthomap GitHub repository <https://github.com/kullrich/orthomap>`_ and at `PyPi] <>`_.

Quick start
-----------

Detailed step-by-step examples are given in the tutorials section.

You can extract an `orthomap` from `OrthoFinder <https://github.com/davidemms/OrthoFinder>`_ results
(for example for zebrafish `Danio rerio` based on the ensembl release-105):

OrthoFinder results for ensembl release-105 can be found here: https://doi.org/10.5281/zenodo.7242264

    >>> from orthomap import of2orthomap, datasets
    >>> # download
    >>> query_orthomap, orthofinder_species_list, of_species_abundance = of2orthomap.get_orthomap(
    ...     seqname='Danio_rerio.GRCz11.cds.longest',
    ...     qt='7955',
    ...     sl='ensembl_105_orthofinder_species_list.tsv',
    ...     oc='ensembl_105_orthofinder_Orthogroups.GeneCount.tsv',
    ...     og='ensembl_105_orthofinder_Orthogroups.tsv',
    ...     continuity=True)

You can show the number of species along query lineage:

    >>> of_species_abundance

You can barplot the counts per taxonomic group (PSname)

    >>> of2orthomap.get_counts_per_ps(query_orthomap).plot.bar(y='counts', x='PSname')

You can get gene to transcript table for `Danio rerio` to be able to match orthomap and scRNA data set gene names:

GTF file here: https://ftp.ensembl.org/pub/release-105/gtf/danio_rerio/Danio_rerio.GRCz11.105.gtf.gz

    >>> query_species_t2g = gtf2t2g.parse_gtf(
    ...     gtf='Danio_rerio.GRCz11.105.gtf.gz',
    ...     g=True, b=True, p=True, v=True, s=True, q=True)

Match

You can calculate transcriptome evolutionary index (TEI) for each cell of a scRNA data set::

   >>> from orthomap import orthomap2tei

You can download/update your local copy of the
NCBI's taxonomy database:

   >>> from orthomap import ncbitax
   >>> ncbitax.update_ncbi()

You can query species lineage information by name:

   >>> from orthomap import qlin
   >>> qlin.get_qlin(q = 'Danio rerio')

or by `taxid`:

   >>> qlin.get_qlin(qt = '7955')

You can also use `orthomap` via the command line::

   $ python src/orthomap/qlin.py -q "Danio rerio"

Table of Contents
-----------------

.. toctree::
   :maxdepth: 2


   tutorials/index.rst
   how-to/index.rst
   references/index.rst

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

.. |GitHub Workflow Status| image:: https://img.shields.io/github/actions/workflow/status/kullrich/orthomap/build_check.yml?branch=main
   :target: https://github.com/kullrich/orthomap/actions/workflows/build_check.yml
.. |PyPI| image:: https://img.shields.io/pypi/v/orthomap?color=blue
   :target: https://pypi.org/project/orthomap/
.. |PyPI - Python Version| image:: https://img.shields.io/pypi/pyversions/orthomap
   :target: https://pypi.org/project/orthomap/
.. |PyPI - Wheel| image:: https://img.shields.io/pypi/wheel/orthomap
   :target: https://pypi.org/project/orthomap/
.. |Licence| image:: https://img.shields.io/badge/License-GPLv3-blue.svg
   :target: https://www.gnu.org/licenses/gpl-3.0
.. |ReadTheDocs| image:: https://readthedocs.org/projects/orthomap/badge/?version=latest
   :target: https://orthomap.readthedocs.io/en/latest/?badge=latest
