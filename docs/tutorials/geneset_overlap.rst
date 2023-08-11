.. _geneset_overlap:

Map OrthoFinder gene names and scRNA gene/transcript names
==========================================================

To be able to link gene ages assignments from an orthomap and gene or transcript
of scRNA dataset, one needs to check the overlap of the annotated gene names.
With the `gtf2t2g` submodule of `orthomap` and the `gtf2t2g.parse_gtf()`
function, one can extract gene and transcript names from a given gene feature
file.