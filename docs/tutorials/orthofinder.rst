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

The species peptide file can be pre-processed e.g. to just contain the longest isoform per gene. (see e.g. `CRBHits::isoform2longest <https://github.com/kullrich/CRBHits>`_.

