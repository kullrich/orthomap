.. _commandline-plaza2orthomap:

Command line function - plaza2orthomap
======================================

::

    orthomap plaza2orthomap -h

    usage: orthomap <sub-command> plaza2orthomap [-h] [-qt QT] [-sl SL] [-og OG] [-out OUT] [-overwrite OVERWRITE]

    optional arguments:
      -h, --help            show this help message and exit
      -qt QT                query species taxID (e.g. use <orthomap qlin -h> to get taxID)
      -sl SL                specify PLAZA species information file <species_information.csv>
      -og OG                specify PLAZA gene family file <genefamily_data.ORTHOFAM.csv> or genefamily_data.HOMFAM.csv
      -out OUT              specify output file <orthomap.tsv> (default: orthomap.tsv)
      -overwrite OVERWRITE  specify if existing output file should be overwritten (default: True)

    plaza2orthomap example:

        # using Orthologous gene family
        $ plaza2orthomap -qt 3702 \
          -sl species_information.csv \
          -og genefamily_data.ORTHOFAM.csv \
          -out 3702.orthofam.orthomap

        # using Homologous gene family
        $ plaza2orthomap -qt 3702 \
          -sl species_information.csv \
          -og genefamily_data.HOMFAM.csv \
          -out 3702.homfam.orthomap