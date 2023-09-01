.. _commandline-of2orthomap:

Command line function - of2orthomap
===================================

::

    orthomap of2orthomap -h

    usage: orthomap <sub-command> of2orthomap [-h] [-seqname SEQNAME] [-qt QT] [-sl SL] [-oc OC] [-og OG] [-out OUT]
                                              [-overwrite OVERWRITE]

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
        $ wget https://zenodo.org/record/7796253/files/ensembl_105_orthofinder_Orthogroups.GeneCount.tsv.zip
        $ wget https://zenodo.org/record/7796253/files/ensembl_105_orthofinder_Orthogroups.tsv.zip
        $ wget https://zenodo.org/record/7796253/files/ensembl_105_orthofinder_species_list.tsv

        # extract orthomap:
        $ of2orthomap -seqname Danio_rerio.GRCz11.cds.longest -qt 7955 \
          -sl ensembl_105_orthofinder_species_list.tsv \
          -oc ensembl_105_orthofinder_Orthogroups.GeneCount.tsv.zip \
          -og ensembl_105_orthofinder_Orthogroups.tsv.zip