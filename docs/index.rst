Welcome to orthomap's documentation!
====================================

`orthomap` is a python package that can:

- extract orthologous maps from OrthoFinder results.

- query species lineage information from the NCBI taxonomic database.

- extract and match gene and transcript information with scRNA data sets.

- calculate and visulaize transcriptome evolutionary index (TEI) on scRNA data sets.

Quick start
-----------

You can extract an orthomap from OrthoFinder results (for example for `Danio rerio` based on the
ensembl release-105)::

   >>> from orthomap import orthomap
   >>> omap = orthomap.get_orthomap(qname = 'Danio_rerio.GRCz11.cds.longest',
   ...        qt = '7955',
   ...        sl = 'ensembl_105_orthofinder_species_list.tsv',
   ...        oc = 'ensembl_105_orthofinder_Orthogroups.GeneCount.tsv',
   ...        og = 'ensembl_105_orthofinder_Orthogroups.tsv')


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
