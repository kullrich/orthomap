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
   **OrthoFinder by default use diamond as the sequence search engine.** To increase sequence search sensitivity, at least use the '-S diamond_ultra_sens' option.
   Even better change the OrthoFinder 'config.josn' and add the option '-k0' for diamond to not only report the best 25, but all sequence search hits.

To change the `'config.json' <https://raw.githubusercontent.com/davidemms/OrthoFinder/master/scripts_of/config.json>`_ and the 'diamond_ultra_sens' option from OrthoFinder, please change the 'cofig.josn' as follows:

   ::

      wget

      "diamond_ultra_sens":{
      "program_type": "search",
      "db_cmd": "diamond makedb --ignore-warnings --in INPUT -d OUTPUT",
      "search_cmd": "diamond blastp --ignore-warnings -k0 -d DATABASE -q INPUT -o OUTPUT --ultra-sensitive -p 1 --quiet -e 0.001 --compress 1"
      },
