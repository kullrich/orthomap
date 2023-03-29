.. _mytai_function:

Correspondance of myTAI and orthomap function
=============================================

::

    library(myTAI)
    data(PhyloExpressionSetExample)
    TAI(PhyloExpressionSetExample)

    >>> import scanpy as sc
    >>> import pandas as pd
    >>> import matplotlib.pyplot as plt
    >>> import seaborn as sns
    >>> from orthomap import of2orthomap, orthomap2tei
    >>> adata = sc.read('PhyloExpressionSetExample.h5ad')
    >>> orthomap2tei.get_tei(
    >>>     adata=adata,
    >>>     gene_id=adata.var.index,
    >>>     gene_age=adata.var['Phylostrata'])

    PlotDistribution(PhyloExpressionSetExample)

    >>> query_orthomap = pd.DataFrame(adata.var.index,
    >>>     columns=['GeneID'])
    >>> query_orthomap['Phylostrata']=adata.var['Phylostrata'].values
    >>> of2orthomap.get_counts_per_ps(omap_df=query_orthomap,
    >>>     psnum_col='Phylostrata',
    >>>     pstaxid_col=None,
    >>>     psname_col=None).plot.bar(y='counts', x='Phylostrata')
    >>> plt.show()

    REMatrix(PhyloExpressionSetExample)
    PlotRE(PhyloExpressionSetExample, Groups=list(1:12))

    >>> rematrix = orthomap2tei.get_rematrix(
    >>>     adata=adata,
    >>>     gene_id=adata.var.index,
    >>>     gene_age=adata.var['Phylostrata'],
    >>>     standard_scale=0)
    >>> rematrix.transpose().plot.line(cmap='Accent')

    pmatrix <- pMatrix(PhyloExpressionSetExample)
    pmatrix
    boxplot(pmatrix, outline=FALSE)

    >>> pmatrix = orthomap2tei.get_pmatrix(
    >>>     adata=adata,
    >>>     gene_id=adata.var.index,
    >>>     gene_age=adata.var['Phylostrata'])
    >>> pd.DataFrame(pmatrix.layers['pmatrix'],
    >>>     index=pmatrix.obs.index).transpose().boxplot(showfliers=False)
    >>> plt.show()

    marker_expression <- PlotGeneSet(ExpressionSet = PhyloExpressionSetExample,
        gene.set = PhyloExpressionSetExample[1:5, 2],
        get.subset = TRUE)
    PlotGeneSet(ExpressionSet = PhyloExpressionSetExample,
        gene.set = PhyloExpressionSetExample[1:5, 2])

    >>> marker_genes = adata.var_names[:5]
    >>> pd.DataFrame(adata[:, marker_genes].X, columns=marker_genes, index=)
