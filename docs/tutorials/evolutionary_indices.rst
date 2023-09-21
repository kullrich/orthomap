.. _evolutionary_indices:

Step 4 - other evolutionary indices
===================================

The phylogenetically based transcriptome evolutionary index (TEI) is one way to weight
transcriptome data. Other evolutionary indices can be used to
weight expression not using a gene age class but other gene based measurements.

The mandatory part for theses indices is a gene based
measurement to be able to assign each gene with a value which is used to weigh its expression.
These values can be discrete or continuous values. Continuous values can be binned first and
used as gene groups to weigh expression.

Other evolutionary indices can be e.g.:

- Tajima'sD
- Nucleotide diversity (within species)
- Nucleotide divergence (between species)
- F-statistics

Here, the TEI measure represents the weighted arithmetic mean
(expression levels as weights for the gene based measurement)
over all categories.

.. math::

    TEI_s = \sum (e_{is} * m_i) / \sum e_{is}

where :math:`TEI_s` denotes the TEI value in developmental stage
:math:`s, e_{is}` denotes the gene expression level of gene :math:`i`
in stage :math:`s`, and :math:`m_i` denotes the corresponding measurement
of gene :math:`i, i = 1,...,N` and :math:`N = total\ number\ of\ genes`.

.. note::
   Please have a look at `Liu & Robinson-Rechavi, 2018 <https://doi.org/10.1093/gbe/evy177>`_ to get more insides how expression transformation can influence TEI calculation.
   By default, the expression values are pre-processed. This includes normalizing the counts per cell (see option `normalize_total`) to a total count of 1e6 (which corresponds to CPM; see option `target_sum`) and log-transformed (see option `log1p`).

How to add other evolutionary indices to scRNA data
---------------------------------------------------

Please download the notebooks from `here <https://raw.githubusercontent.com/kullrich/orthomap/main/docs/notebooks/evolutionary_indices.ipynb>`_
or please click below to view the content.

.. toctree::

   ../notebooks/evolutionary_indices