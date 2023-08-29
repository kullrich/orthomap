.. _tutorial:

Tutorials
=========

This section contains a variety of tutorials that should help get you started
with the `orthomap` package.

What the tutorial covers
------------------------

orthomap - steps
----------------

- overview picture goes here

This section contains the main steps of `orthomap` to extract gene age information for a query species up to linking
the extracted gene age classes and expression data of single-cell data sets.

- :doc:`orthofinder`: This tutorial introduces how to run your own OrthoFinder analysis (step 0).
- :doc:`query_lineage`: This tutorial introduces how to get taxonomic information (step 1).
- :doc:`get_orthomap`: This tutorial introduces how to extract an orthomap (gene age class) from OrthoFinder results (step 2).
- :doc:`geneset_overlap`: This tutorial introduces how to match gene or transcript IDs between an orthomap and scRNA data (step 3).
- :doc:`add_tei`: This tutorial introduces how to add a transcriptome evolutionary index (short: TEI) to scRNA data (step 4).

orthomap - downstream analysis (step 5)
---------------------------------------

- :doc:`pstrata`: This tutorial introduces partial TEI and its contribution to the global TEI per cell or cell type (step 5).
- :doc:`plotting`: This tutorial introduces some basic concepts of plotting results (step 5).
- :doc:`relative_expression`: This tutorial introduces relative expression per gene age class and its contribution to the global TEI per cell or cell type (step 5).

case studies
------------

- :doc:`nematode_example`: This notebook covers a re-analysis of nematode (Caenorhabditis elegans) scRNA data.
- :doc:`zebrafish_example`: This notebook covers a re-analysis of zebrafish (Danio rerio) scRNA data.
- :doc:`frog_example`: This notebook covers a re-analysis of frog (Xenopus tropicalis) scRNA data.
- :doc:`mouse_example`: This notebook covers a re-analysis of mouse (Mus musculus) scRNA data.
- :doc:`hydra_example`: This notebook covers a re-analysis of hydra (Hydra vulgaris) scRNA data.

command

- :doc:`commandline`: This tutorial covers which orthomap functions can be run via the command line.
- :doc:`mytai`: This tutorial covers which orthomap functions correspond to myTAI functions.

.. note::
   A demo dataset is available for each of the tutorial notebooks above.
   These datasets allow you to begin exploring `orthomap` even if you do not have any data at any step in the analysis
   pipeline.

Prerequisites
-------------

- This tutorial assumes that you have basic **Python programming experience**.
  In particular, we assume you are familiar with using a notebook from the following python data science libraries:
  **jupyter**.
- To better understand plotting and data access, the user should try to get familiar with the python libraries:
  **pandas**, **matplotlib** and **seaborn**.
- `orthomap` is a python package but part of it can be run on the command line. For the installation of `orthomap`,
  we recommend using `Anaconda <https://anaconda.org>`_
  (`see here <https://orthomap.readthedocs.io/en/latest/installation/index.html>`_).
  If you are not familiar with Anaconda or python environment management,
  please use `our pre-built docker image <https://orthomap.readthedocs.io/en/latest/installation/index.html#docker-image>`_.

Code and data availability
--------------------------

- We provide links for the notebook in each section.

- You can download the demo input data using `orthomap` data loading function in the notebooks.
  `see here <https://orthomap.readthedocs.io/en/latest/modules/orthomap.html#modules-for-dataset-downloads>`_)

Getting started
---------------

If you are running `orthomap` for the first time, we recommend getting started with the :doc:`nematode_example` and
then, proceed to :doc:`zebrafish_example`.

Table of Contents
-----------------

.. toctree::
   :caption: orthomap - steps
   :maxdepth: 1

   orthofinder
   query_lineage
   get_orthomap
   geneset_overlap
   add_tei

.. toctree::
   :caption: Command line
   :maxdepth: 1

   commandline.cds2aa
   commandline.gtf2t2g
   commandline.ncbitax
   commandline.of2orthomap
   commandline.plaza2orthomap
   commandline.qlin


.. toctree::
   :caption: myTAI - function correspondance
   :maxdepth: 1

   mytai
