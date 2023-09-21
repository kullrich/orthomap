.. _add_tei:

Step 4 - TEI calculation
========================

The phylogenetically based transcriptome evolutionary index (TEI) is similar to
`Domazet-Loso & Tautz, 2010 <https://doi.org/10.1038/nature09632>`_.

The TEI measure represents the weighted arithmetic mean
(expression levels as weights for the phylostratum value)
over all evolutionary age categories denoted as *phylostra*.

.. math::

    TEI_s = \sum (e_{is} * ps_i) / \sum e_{is}

where :math:`TEI_s` denotes the TEI value in developmental stage
:math:`s, e_{is}` denotes the gene expression level of gene :math:`i`
in stage :math:`s`, and :math:`ps_i` denotes the corresponding phylostratum
of gene :math:`i, i = 1,...,N` and :math:`N = total\ number\ of\ genes`.

.. note::
   Please have a look at `Liu & Robinson-Rechavi, 2018 <https://doi.org/10.1093/gbe/evy177>`_ to get more insides how expression transformation can influence TEI calculation.
   By default, the expression values are pre-processed. This includes normalizing the counts per cell (see option `normalize_total`) to a total count of 1e6 (which corresponds to CPM; see option `target_sum`) and log-transformed (see option `log1p`).

How to add a transcriptome evolutionary index (short: TEI) to scRNA data
------------------------------------------------------------------------

Please download the notebooks from `here <https://raw.githubusercontent.com/kullrich/orthomap/main/docs/notebooks/add_tei.ipynb>`_
or please click below to view the content.

.. toctree::

   ../notebooks/add_tei