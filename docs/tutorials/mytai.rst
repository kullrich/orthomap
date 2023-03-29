.. _mytai_function:

Correspondance of myTAI and orthomap function
=============================================

`myTAI::TAI()` in R

.. code-block:: R

    install.packages("myTAI")
    library(myTAI)
    data(PhyloExpressionSetExample)
    TAI(PhyloExpressionSetExample)

`orthomap2tei.get_tei()` in Python

.. code-block:: Python

    import scanpy as sc
    import pandas as pd
    import matplotlib.pyplot as plt
    import seaborn as sns
    from orthomap import of2orthomap, orthomap2tei, datasets
    adata = datasets.mytai_example(datapath='.')
    orthomap2tei.get_tei(
        adata=adata,
        gene_id=adata.var.index,
        gene_age=adata.var['Phylostrata'])

`myTAI::PlotDistribution()` in R

.. code-block:: R

    PlotDistribution(PhyloExpressionSetExample)

`of2orthomap.get_counts_per_ps()` in Python

.. code-block:: Python

    query_orthomap = pd.DataFrame(adata.var.index,
        columns=['GeneID'])
    query_orthomap['Phylostrata']=adata.var['Phylostrata'].values
    of2orthomap.get_counts_per_ps(omap_df=query_orthomap,
        psnum_col='Phylostrata',
        pstaxid_col=None,
        psname_col=None).plot.bar(y='counts', x='Phylostrata')
    plt.show()

`myTAI::REMatrix()` and `myTAI::PlotRE()` in R

.. code-block:: R

    REMatrix(PhyloExpressionSetExample)
    PlotRE(PhyloExpressionSetExample, Groups=list(1:12))

`orthomap2tei.get_rematrix()` in Python

.. code-block:: Python

    rematrix = orthomap2tei.get_rematrix(
        adata=adata,
        gene_id=adata.var.index,
        gene_age=adata.var['Phylostrata'],
        standard_scale=0)
    rematrix.transpose().plot.line(cmap='Accent')
    plt.show()

`myTAI::pStrata()` and `myTAI::PlotContribution()` in R

.. code-block:: R

    pstrata <- pStrata(PhyloExpressionSetExample)
    PlotContribution(PhyloExpressionSetExample, "PS")

`orthomap2tei.get_pstrata()` in Python

.. code-block:: Python

    pstrata = orthomap2tei.get_pstrata(
        adata=adata,
        gene_id=adata.var.index,
        gene_age=adata.var['Phylostrata'])
    pstrata[0]
    pstrata[0].transpose().plot.line(cmap='Accent', stacked=True)
    plt.show()

`myTAI::pMatrix()` in R

.. code-block:: R

    pmatrix <- pMatrix(PhyloExpressionSetExample)
    pmatrix
    boxplot(pmatrix, outline=FALSE)

`orthomap2tei.get_pmatrix()` in Python

.. code-block:: Python

    pmatrix = orthomap2tei.get_pmatrix(
        adata=adata,
        gene_id=adata.var.index,
        gene_age=adata.var['Phylostrata'])
    pd.DataFrame(pmatrix.layers['pmatrix'].toarray(),
        index=pmatrix.obs.index).transpose().boxplot(showfliers=False)
    plt.show()

`myTAI::PlotGeneSet()` in R

.. code-block:: R

    marker_expression <- PlotGeneSet(ExpressionSet = PhyloExpressionSetExample,
        gene.set = PhyloExpressionSetExample[1:5, 2],
        get.subset = TRUE)
    PlotGeneSet(ExpressionSet = PhyloExpressionSetExample,
        gene.set = PhyloExpressionSetExample[1:5, 2])

`scanpy` in Python

.. code-block:: Python

    marker_genes = adata.var_names[:5]
    marker_expression = pd.DataFrame(adata[:, marker_genes].X.toarray(),
        columns=marker_genes, index=adata.obs.index)
    marker_expression.plot.line(cmap='Accent')
    plt.show()

`myTAI::PlotMeans()` in R

.. code-block:: R

    PlotMeans(PhyloExpressionSetExample, Groups=list(1:12))

`orthomap2tei.get_ematrix()` in Python

.. code-block:: Python

    ematrix = orthomap2tei.get_ematrix(
        adata=adata,
        group_by_var='Phylostrata')
    ematrix.transpose().plot.line(cmap='Accent')
    plt.show()
