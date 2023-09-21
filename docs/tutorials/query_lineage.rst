.. _query_lineage:

Step 1 - get taxonomic information
==================================

Given a species name or taxonomic ID, the query species lineage information is
extracted with the help of the `ete3 python toolkit <http://etetoolkit.org/>`_
and the `NCBI taxonomy <https://www.ncbi.nlm.nih.gov/taxonomy>`_ (Huerta-Cepas et al., 2016).

This information is needed alongside with the taxonomic classifications for all species used in the OrthoFinder comparison.

.. note::
   The NCBI taxonomy database needs to be downloaded using the package `ete3`.
   A parsed version of it will be stored at your home directory: `~/.etetoolkit/taxa.sqlite`.

`orthomap` provides a command line function to download or upgrade the NCBI taxonomy database (see :ref:`commandline.ncbitax`).

How to get taxonomic information about your query species
---------------------------------------------------------

.. note::
   The orthomap submodule `qlin <https://orthomap.readthedocs.io/en/latest/modules/orthomap.html#modules-for-query-lineage>`_
   helps to get taxonomic information for you with the `qlin.get_qlin() <https://orthomap.readthedocs.io/en/latest/modules/orthomap.qlin.html#orthomap.qlin.get_qlin>`_ function.
   The function will output the following list:
   - query name `query_lineage[0]`
   - query taxID `query_lineage[1]`
   - query lineage `query_lineage[2]`
   - query lineage dictionary `query_lineage[3]`
   - query lineage zip `query_lineage[4]`
   - query lineage names `query_lineage[5]`
   - reverse query lineage `query_lineage[6]`
   - query kingdom `query_lineage[7]`


Please download the notebooks from `here <https://raw.githubusercontent.com/kullrich/orthomap/main/docs/notebooks/query_lineage.ipynb>`_
or please click below to view the content.

.. toctree::

   ../notebooks/query_lineage