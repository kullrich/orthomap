.. _orthofinder:

Step 0 - run OrthoFinder
========================

In order to extract an orthomap from `OrthoFinder <https:https://github.com/davidemms/OrthoFinder>`_ results, one needs to run `OrthoFinder <https:https://github.com/davidemms/OrthoFinder>`_.

Mandatory OrthoFinder results files
-----------------------------------

To be able to extract gene age classes from `OrthoFinder <https:https://github.com/davidemms/OrthoFinder>`_ results one needs the following two files:

- <Orthogroups.GeneCount.tsv> or <Orthogroups.GeneCount.tsv.zip>
- <Orthogroups.tsv> or <Orthogroups.tsv.zip>

.. note::
   It is not necessary to run the full analysis of OrthoFinder. Since only the <Orthogroups.GeneCount.tsv> and <Orthogroups.tsv>
   files are needed, one can start OrthoFinder with the `-og` option and skip further analysis steps.

Mandatory species information
-----------------------------

To be able to extract gene age classes for a given query species the sequence names used for the `OrthoFinder <https:https://github.com/davidemms/OrthoFinder>`_
analysis needs to be listed in a tab delimited file (called species list) with the following two columns:

- <OrthoFinder name>
- <species taxID>

One line of the species list file should look like this:

<OrthoFinder name><tab><species taxID>

In total all species from the `OrthoFinder <https:https://github.com/davidemms/OrthoFinder>`_ run that should be taken
into account for the gene age class assignment need to be listed in that file.

.. note::
   The <OrthoFinder name> is not the common species name but the sequence (FASTA file) name used to start OrthoFinder.

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

To extract the longest isoform `orthomap`

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

