.. _paul15_example:

hematopoiesis - case study
==========================

Reconstructing myeloid and erythroid differentiation using mouse scRNA data of `Paul et al. (2015) <http://doi.org/10.1016/j.cell.2015.11.013>`_.

see also:

- https://scanpy-tutorials.readthedocs.io/en/latest/paga-paul15.html
- https://morris-lab.github.io/CellOracle.documentation/notebooks/03_scRNA-seq_data_preprocessing/scanpy_preprocessing_with_Paul_etal_2015_data.html

*Mus musculus* hematopoiesis single-cell data analysis example
--------------------------------------------------------------

.. note::
   Since in the scRNA data only gene names are provided and not gene/transcript IDs, there is only partial overlap
   with the general transfer format (GTF) from mouse (`Ensembl - release 110 - Mus_musculus.GRCm39.110.gtf.gz <https://ftp.ensembl.org/pub/release-110/gtf/mus_musculus/Mus_musculus.GRCm39.110.gtf.gz>`_)
   and as such also with the OrthoFinder results (3044/3451). This is why it should be preferred to at least provide an additional column in the adata.var with a common GeneID from either NCBI or ensembl.


Please download the notebooks from `here <https://raw.githubusercontent.com/kullrich/orthomap/main/docs/notebooks/paul15_example.ipynb>`_ .
Or please click below to view the content.

.. toctree::

   ../notebooks/paul15_example