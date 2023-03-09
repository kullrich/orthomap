orthomap command line tools
===========================

orthomap
--------

::

    orthomap -h

    options:
      -h, --help            show this help message and exit

    sub-commands:
      {gtf2t2g,ncbitax,of2orthomap,qlin}
                            sub-commands help
        gtf2t2g             extracts transcript to gene table from GTF <gtf2t2g -h>
        ncbitax             update local ncbi taxonomy database <ncbitax -h>
        of2orthomap         extract orthomap from OrthoFinder output for query species<orthomap -h>
        qlin                get query lineage based on ncbi taxonomy <qlin -h>

orthomap gtf2t2g
^^^^^^^^^^^^^^^^

::

    orthomap gtf2t2g -h

    usage: orthomap <sub-command> gtf2t2g [-h] [-i I] [-o O] [-g] [-b] [-p] [-v] [-s] [-q] [-overwrite OVERWRITE]

    options:
      -h, --help            show this help message and exit
      -i I                  specify GTF input file
      -o O                  specify output file [optional]
      -g                    specify if gene names should be appended if they exist
      -b                    specify if gene biotype should be appended if they exist
      -p                    specify if protein id should be appended if they exist
      -v                    specify if gene/transcript/protein version should be appended
      -s                    specify if summary should be printed
      -q                    specify if output should be quite
      -overwrite OVERWRITE  specify if existing output file should be overwritten (default: True)

    gtf2t2g example:

        # to get GTF from Mus musculus on Linux run:
        $ wget https://ftp.ensembl.org/pub/release-108/gtf/mus_musculus/Mus_musculus.GRCm39.108.chr.gtf.gz

        # on Mac:
        $ curl https://ftp.ensembl.org/pub/release-108/gtf/mus_musculus/Mus_musculus.GRCm39.108.chr.gtf.gz --remote-name

        # create t2g from GTF:
        $ gtf2t2g -i Mus_musculus.GRCm39.108.chr.gtf.gz -o Mus_musculus.GRCm39.108.chr.gtf.t2g.tsv -g -b -p -v -s

orthomap ncbitax
^^^^^^^^^^^^^^^^

::

    orthomap ncbitax -h

    usage: orthomap <sub-command> ncbitax [-h] [-u]

    options:
      -h, --help  show this help message and exit
      -u          update

    ncbitax example:

        #update ncbi taxonomy database:
        ncbitax -u

orthomap of2orthomap
^^^^^^^^^^^^^^^^^^^^

::

    orthomap of2orthomap -h

    usage: orthomap <sub-command> of2orthomap [-h] [-seqname SEQNAME] [-qt QT] [-sl SL] [-oc OC] [-og OG] [-out OUT] [-overwrite OVERWRITE]

    options:
      -h, --help            show this help message and exit
      -seqname SEQNAME      sequence name of the query species in orthofinder(see column names of <Orthogroups.tsv>)
      -qt QT                query species taxid (e.g. use <orthomap qlin -h> to get taxid)
      -sl SL                species list as <orthofinder name><tab><species taxid> (only samples in this list will be processed)
      -oc OC                specify orthofinder <Orthogroups.GeneCounts.tsv> (see Orthogroups directory)
      -og OG                specify orthofinder <Orthogroups.tsv> (see Orthogroups directory)
      -out OUT              specify output file <orthomap.tsv> (default: orthomap.tsv)
      -overwrite OVERWRITE  specify if existing output file should be overwritten (default: True)

    of2orthomap example:

        # download OrthoFinder example:
        $ wget https://github.com/kullrich/orthomap/raw/main/examples/ensembl_105_orthofinder_Orthogroups.GeneCount.tsv.zip
        $ wget https://github.com/kullrich/orthomap/raw/main/examples/ensembl_105_orthofinder_Orthogroups.tsv.zip
        $ wget https://github.com/kullrich/orthomap/raw/main/examples/ensembl_105_orthofinder_species_list.tsv

        # extract orthomap:
        $ of2orthomap -seqname Danio_rerio.GRCz11.cds.longest -qt 7955 \
          -sl ensembl_105_orthofinder_species_list.tsv \
          -oc ensembl_105_orthofinder_Orthogroups.GeneCount.tsv.zip \
          -og ensembl_105_orthofinder_Orthogroups.tsv.zip

orthomap qlin
^^^^^^^^^^^^^

::

    orthomap qlin -h

    usage: orthomap <sub-command> qlin [-h] [-q Q] [-qt QT]

    options:
      -h, --help  show this help message and exit
      -q Q        query species name
      -qt QT      query species taxid

    qlin example:

        # get query lineage to be used with orthomap later on using query species taxid
        # Mus musculus; 10090
        $ qlin -qt 10090

        # using query species name
        $ qlin -q "Mus musculus"

Modules for dataset downloads
=============================

 .. toctree::

    orthomap.datasets

Modules for eggnog
==================

 .. toctree::

    orthomap.eggnog2orthomap

Modules for GTF handling
========================

 .. toctree::

    orthomap.gtf2t2g

Modules for NCBI taxonomy
=========================

 .. toctree::

    orthomap.ncbitax

Modules for OrthoFinder
=======================

 .. toctree::

    orthomap.of2orthomap

Modules for single-cell data
============================

 .. toctree::

    orthomap.orthomap2tei

Modules for query lineage
=========================

 .. toctree::

    orthomap.qlin
