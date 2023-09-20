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
   You can see the sequence names used by OrthoFinder in the header of either the <Orthogroups.GeneCount.tsv> or <Orthogroups.tsv> file.

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

.. note::
   If you have used a conda environment to install `OrthoFinder <https:https://github.com/davidemms/OrthoFinder>`_,
   the 'config.json' of `OrthoFinder <https:https://github.com/davidemms/OrthoFinder>`_
   can be found in in the location '~./conda/envs/ENVNAME/bin/scripts_of/config.json'.

To change the `'config.json' <https://raw.githubusercontent.com/davidemms/OrthoFinder/master/scripts_of/config.json>`_ and the 'diamond_ultra_sens' option from OrthoFinder, please change the 'cofig.json' as follows:

   ::

      "diamond_ultra_sens":{
      "program_type": "search",
      "db_cmd": "diamond makedb --ignore-warnings --in INPUT -d OUTPUT",
      "search_cmd": "diamond blastp --ignore-warnings -k0 -d DATABASE -q INPUT -o OUTPUT --ultra-sensitive -p 1 --quiet -e 0.001 --compress 1"
      },


Use LAST with OrthoFinder
-------------------------

To use e.g.: `last <https://gitlab.com/mcfrith/last>`_ as the sequence search engine, please install it and make the executables available to your PATH.

The source code can be either accessed here:

- `last source code <https://gitlab.com/mcfrith/last/-/tags>`_

or you might want to install with bioconda:

- `bioconda last <https://anaconda.org/bioconda/last>`_

To use `last <https://gitlab.com/mcfrith/last>`_ as a new sequence serach engine,
please change the 'config.json' as follows:

   ::

      "last":{
      "program_type": "search",
      "db_cmd": "lastdb -p -cR01 OUTPUT INPUT",
      "search_cmd": "lastal -f BlastTab+ -D 1e6 DATABASE INPUT | sed -n '/^#/!p' > OUTPUT"
      },


Typical run command
-------------------

- using diamond

   ::

      orthofinder -t 32 -a 8 -og -o diamond_output/ -S diamond_ultra_sens -f folder_with_peptides/


- using last

   ::

      orthofinder -t 32 -a 8 -og -o last_output/ -S last -f folder_with_peptides/


Adding a new species to an existing OrthoFinder result
------------------------------------------------------

Working with pre-annotated scRNA data is sometimes cumbersome, since an "older" genome annotation version
was used for your species of interest. It might be that an original published study used transcriptome information and
not genome annotation.

In both cases it might be difficult to find 1-to-1 sequence hits between the "older" and "newer" gene annotation version.

However, with `OrthoFinder <https:https://github.com/davidemms/OrthoFinder>`_ it is possible to add a new species to an existing analysis
which can be helpful in this situation.

Here, a short proposal is given how to deal with that situation.
In the original publication of `Plass, Solana et al, 2018 <https://doi.org/10.1126/science.aaq1723>`_
and the planaria species *Schmidtea mediterranea* `GSE103633<https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103633>`_
a transcriptome was used. However, the exists an annotated genome for the same species
(`Schmidtea mediterranea PRJNA12585 peptides <https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS18/species/schmidtea_mediterranea/PRJNA12585/schmidtea_mediterranea.PRJNA12585.WBPS18.protein.fa.gz>`_), but the scRNA data uses
the transcriptome contig IDs as gene IDs.

**Proposed workflow:**

- call ORFs/CDS from the given transcriptome

ORF/CDS extraction can be done with e.g. `TransDecoder<https://github.com/TransDecoder/TransDecoder>`_ or
using `miniprot<https://github.com/lh3/miniprot>`_ with the "newer" annotated peptides followed by `miniprothint<https://github.com/tomasbruna/miniprothint>`_ or
using `GALBA<https://github.com/Gaius-Augustus/GALBA>`_

   ::

       miniprot dd_Smed_v6.pcf.contigs.fasta schmidtea_mediterranea.PRJNA12585.WBPS18.protein.fa --aln > miniprot.aln
       miniprot_boundary_scorer -o miniprot_parsed.gff -s blosum62.csv < miniprot.aln
       miniprothint.py miniprot_parsed.gff --workdir miniprothint


- extract and convert CDS into peptides from the given transcriptome

extraction and direct conversion into peptides can be done with e.g. `gffread<https://github.com/gpertea/gffread>`_

Here, first the original contig IDs are added to the gene IDs so that later a mapping against the scRNA data is possible.

   ::

       awk -F '\t' -vOFS='\t' '{if($3=="mRNA"){gsub("ID=","ID="$1"::",$9)}; if($3!="mRNA"){gsub("Parent=", "Parent="$1"::", $9)}; print $0}' miniprot_parsed.gff > miniprot_parsed_IDs.gff
           gffread -x dd_Smed_v6_miniprot_parsed.x.fasta -y dd_Smed_v6_miniprot_parsed.pep.fasta -g dd_Smed_v6.pcf.contigs.fasta miniprot_parsed_IDs.gff

Now one can use the extracted peptides with `OrthoFinder <https:https://github.com/davidemms/OrthoFinder>`_ to add them to an existing `OrthoFinder <https:https://github.com/davidemms/OrthoFinder>`_ run.

- Place the new species peptide files in a separate folder

   ::

       orthofinder -t 32 -a 8 -og -S last -b last_output/Results_Sep13/WorkingDirectory/ -f new_species/