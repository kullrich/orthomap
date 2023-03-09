orthomap command line tools
===========================

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
