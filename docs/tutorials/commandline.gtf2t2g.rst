.. _gtf2t2g_cmd:

Command line - gtf2t2g
======================

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
