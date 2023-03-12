.. _orthofinder:

How to run OrthoFinder
======================

In order to extract an orthomap from OrthoFinder results, one needs to run `OrthoFinder <https:https://github.com/davidemms/OrthoFinder>`_.

Install OrthoFinder
-------------------

OrthoFinder installation using conda
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

  ::

      conda install -c bioconda orthofinder

Run OrthoFinder
---------------

OrthoFinder can take a folder as input which should contain all input fasta files, one file per species.

The species peptide file can be pre-processed e.g. to just contain the longest isoform per gene.

  ::

      wget https://ftp.ensembl.org/pub/release-105/fasta/danio_rerio/cds/Danio_rerio.GRCz11.cds.all.fa.gz
      gunzip Danio_rerio.GRCz11.cds.all.fa.gz
      orthomap cds2aa -i Danio_rerio.GRCz11.cds.all.fa -r ENSEMBL -o Danio_rerio.GRCz11.aa.all.longest.fa

.. warning::
   **OrthoFinder by default use diamond as the sequence search engine.** Please to increase sequence search sensitivity, at least use the '-S diamond_ultra_sens' option.
   Even better change the OrthoFinder 'config' and add the option '-k0' for diamond to not only report the best 25, but all search hits.