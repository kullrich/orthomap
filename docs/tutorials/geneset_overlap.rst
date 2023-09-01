.. _geneset-overlap:

Step 3 - map gene/transcript IDs
================================

To be able to link gene ages assignments from an orthomap and gene or transcript
of scRNA dataset, one needs to check the overlap of the annotated gene names.

With the `gtf2t2g` submodule of `orthomap` and the `gtf2t2g.parse_gtf()`
function, one can extract gene and transcript names from a given gene feature
file.

How to match gene or transcript IDs between an orthomap and scRNA data
----------------------------------------------------------------------
