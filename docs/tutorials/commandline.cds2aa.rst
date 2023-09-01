.. _commandline-cds2aa:

Command line function - cds2aa
==============================

::

    orthomap cds2aa -h

    usage: orthomap <sub-command> cds2aa [-h] [-i I] [-o O] [-t T] [-r R]

    options:
      -h, --help  show this help message and exit
      -i I        specify fasta input file
      -o O        specify output file [optional]
      -t T        transtable [default: std]
      -r R        specify CDS source to retain longest isoform

    cds2aa example:

        # to get CDS from Danio rerio on Linux run:
        $ wget https://ftp.ensembl.org/pub/release-105/fasta/danio_rerio/cds/Danio_rerio.GRCz11.cds.all.fa.gz
        $ gunzip Danio_rerio.GRCz11.cds.all.fa.gz

        # on Mac:
        $ curl https://ftp.ensembl.org/pub/release-105/fasta/danio_rerio/cds/Danio_rerio.GRCz11.cds.all.fa.gz --remote-name
        $ gunzip Danio_rerio.GRCz11.cds.all.fa.gz

        # translate and retain longest isoform from CDS fasta file:
        $ cds2aa -i Danio_rerio.GRCz11.cds.all.fa -r ENSEMBL -o Danio_rerio.GRCz11.aa.all.longest.fa