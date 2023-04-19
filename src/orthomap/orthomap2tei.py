"""
Author: Kristian K Ullrich
date: April 2023
email: ullrich@evolbio.mpg.de
License: GPL-3
"""


import os
import scipy
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import seaborn as sns
from alive_progress import alive_bar
from statannot import add_stat_annotation


def read_orthomap(orthomapfile):
    """
    This function reads a pre-calculated orthomap file <GeneID><tab><Phylostratum>.

    :param orthomapfile: File name of pre-calculated orthomap file.
    :return: Orthomap.

    :type orthomapfile: str
    :rtype: pandas.DataFrame

    Example
    -------
    >>> from orthomap import orthomap2tei, datasets
    >>> # download pre-calculated orthomap
    >>> #query_orthomap = orthomap2tei.read_orthomap(orthomapfile='Sun2021_Orthomap.tsv')
    >>> sun21_orthomap_file = datasets.sun21_orthomap(datapath='.')
    >>> # load query species orthomap
    >>> query_orthomap = orthomap2tei.read_orthomap(orthomapfile=sun21_orthomap_file)
    >>> query_orthomap
    """
    orthomap = None
    if os.path.exists(orthomapfile):
        orthomap = pd.read_csv(orthomapfile,
                               delimiter='\t')
    return orthomap


def geneset_overlap(geneset1,
                    geneset2):
    """
    This function shows the overlap of two lists. To check e.g. <GeneID> from an orthomap and <adata.var_names>
    from an AnnData object.

    :param geneset1: List of gene or transcript names set 1.
    :param geneset2: List of gene or transcript names set 2.
    :return: Overlap.

    :type geneset1: list
    :type geneset2: list
    :rtype: pandas.DataFrame

    Example
    -------
    >>> from orthomap import orthomap2tei
    >>> geneset1 = ['g1.1', 'g1.2', 'g2.1', 'g3.1', 'g3.2']
    >>> geneset2 = ['g1.1', 'g2.1', 'g3.1']
    >>> orthomap2tei.geneset_overlap(geneset1, geneset2)
    """
    g1 = set(geneset1)
    g2 = set(geneset2)
    g1_g2 = g1.intersection(g2)
    d = {'g1_g2_overlap': [len(g1_g2)],
         'g1_ratio': [len(g1_g2)/len(g1)],
         'g2_ratio': [len(g1_g2)/len(g2)]}
    df = pd.DataFrame(data=d)
    return df


def replace_by(x_orig,
               xmatch,
               xreplace):
    """
    This function assumes that <x_orig> and <xmatch> match and will return <xreplace> sorted by <x_orig>.
    It is mandatory that <xmatch> and <xreplace> have the same length and reflect pairs:
    <xmatch[0]> is the original value and <xreplace[0]> is the corresponding new value.

    :param x_orig: List of original values to be used for sorting.
    :param xmatch: List of matches to the original values. Each xmatch position pairs with xreplace position.
    :param xreplace: List of replace values. Each xreplace position pairs with xmatch position.
    :return: Replacement ordered by the original values.

    :type x_orig: list
    :type xmatch: list
    :type xreplace: list
    :rtype: list

    Example
    -------
    >>> from orthomap import orthomap2tei
    >>> geneset1 = ['g1.1', 'g1.2', 'g2.1', 'g3.1', 'g3.2']
    >>> geneset2 = ['g1.1', 'g2.1', 'g3.1', 'g5.1']
    >>> transcriptset2 = ['t1.1', 't2.1', 't3.1', 't5.1']
    >>> orthomap2tei.replace_by(
    >>>     x_orig=geneset1,
    >>>     xmatch=geneset2,
    >>>     xreplace=transcriptset2)
    """
    replace_dict = {}
    for i, j in enumerate(xmatch):
        replace_dict[j] = xreplace[i]
    x_new = [replace_dict[x] if x in replace_dict else np.nan for x in x_orig]
    return x_new


def _keep_min_max(df,
                  keep='min',
                  dup_col='GeneID',
                  sort_col='Phylostrata'):
    """
    A helper function to keep either the minimal or maximal value based on a duplication column.

    :param df: Expects orthomap DataFrame, but can be any if dup_col and sort_col are present.
    :param keep: Either define 'min' (ascending pre-sorting) or 'max' (non-ascending pre-sorting) to keep duplicates.
    :param dup_col: Column name(s) to be searched for duplicates.
    :param sort_col: Column names(s) to sort DataFrame.
    :return: Reduced DataFrame.

    :type df: pandas.DataFrame
    :type keep: str
    :type dup_col: str
    :type sort_col: str
    :rtype: pandas.DataFrame

    Example
    -------
    >>> import pandas as pd
    >>> from orthomap import orthomap2tei
    >>> # create artificial DataFrame
    >>> my_orthomap = pd.DataFrame.from_dict(
    >>>     {'GeneID':['g1', 'g1', 'g2', 'g3', 'g3'],
    >>>     'Phylostrata':[3, 1, 2, 5, 7]})
    >>> # keep min value
    >>> orthomap2tei._keep_min_max(
    >>>     df=my_orthomap,
    >>>     keep='min')
    >>> # keep max value
    >>> orthomap2tei._keep_min_max(
    >>>     df=my_orthomap,
    >>>     keep='max')
    """
    df_sorted = None
    if keep == 'min':
        df_sorted = df.sort_values(by=[sort_col],
                                   ascending=True)
    if keep == 'max':
        df_sorted = df.sort_values(by=[sort_col],
                                   ascending=False)
    keep_idx = ~df_sorted.duplicated([dup_col])
    df_out = df_sorted[keep_idx]
    return df_out


def _split_gene_id_by_gene_age(gene_id,
                               gene_age,
                               keep='min',
                               adata=None):
    """
    A helper function to group <GeneID> by <Phylostrata>.

    :param gene_id: Expects GeneID column from orthomap DataFrame.
    :param gene_age: Expects Phylostratum column from orthomap DataFrame.
    :param keep: In case of duplicated GeneIDs with different Phylostrata assignments, either keep 'min' or 'max' value.
    :param adata: AnnData object of shape n_obs × n_vars. Rows correspond to cells and columns to genes.
    :return: Dictonary with gene_age as keys and gene_id as values.

    :type gene_id: list
    :type gene_age: list
    :type keep: str
    :type adata: AnnData
    :rtype: dictionary

    Example
    -------
    >>> from orthomap import orthomap2tei, datasets
    >>> # download pre-calculated orthomap
    >>> #query_orthomap = orthomap2tei.read_orthomap(orthomapfile='Sun2021_Orthomap.tsv')
    >>> sun21_orthomap_file = datasets.sun21_orthomap(datapath='.')
    >>> # load query species orthomap
    >>> query_orthomap = orthomap2tei.read_orthomap(orthomapfile=sun21_orthomap_file)
    >>> orthomap2tei._split_gene_id_by_gene_age(
    >>>     gene_id=query_orthomap['GeneID'],
    >>>     gene_age=query_orthomap['Phylostrata'],
    >>>     keep='min')
    """
    id_age_df = pd.DataFrame(data={'GeneID': gene_id,
                                   'Phylostrata': gene_age})
    # check and drop duplicated GeneID
    id_age_df_keep = _keep_min_max(df=id_age_df,
                                   keep=keep,
                                   dup_col='GeneID',
                                   sort_col='Phylostrata')
    # get overlap
    if adata is not None:
        gene_intersection = pd.Index(id_age_df_keep['GeneID']).intersection(adata.var_names)
        id_age_df_keep_subset = id_age_df_keep.loc[id_age_df_keep['GeneID'].isin(gene_intersection)]
        id_age_df_keep_subset = id_age_df_keep_subset.sort_values('GeneID')
    else:
        id_age_df_keep_subset = id_age_df_keep
    phylostrata = list(set(id_age_df_keep_subset['Phylostrata']))
    gene_id_gene_age_dict = {}
    for pk_idx, pk in enumerate(phylostrata):
        gene_id_gene_age_dict[str(pk)] = list(
            id_age_df_keep_subset[id_age_df_keep_subset['Phylostrata'] == pk]['GeneID'])
    return gene_id_gene_age_dict


def _get_counts(adata,
                layer=None,
                normalize_total=False,
                log1p=False,
                target_sum=1e6):
    """
    A helper function to pre-process AnnData counts.

    :param adata: AnnData object of shape n_obs × n_vars. Rows correspond to cells and columns to genes.
    :param layer: Layer to work on instead of X. If None, X is used.
    :param normalize_total: Normalize counts per cell prior TEI calculation.
    :param log1p: Logarithmize the data matrix prior TEI calculation.
    :param target_sum: After normalization, each observation (cell) has a total count equal to target_sum.
    :return: Processed AnnData counts.

    :type adata: AnnData
    :type layer: str
    :type normalize_total: bool
    :type log1p: bool
    :type target_sum: float
    :rtype: scipy.sparse._csr.csr_matrix

    Example
    -------
    >>> import scanpy as sc
    >>> from orthomap import orthomap2tei, datasets
    >>> # download and load scRNA data
    >>> #packer19_small = sc.read('packer19_small.h5ad')
    >>> packer19_small = datasets.packer19_small(datapath='.')
    >>> orthomap2tei._get_counts(adata=packer19_small)
    """
    adata_counts = adata.X
    if layer is not None:
        adata_counts = adata.layers[layer]
    if normalize_total and log1p:
        adata_norm = sc.pp.normalize_total(adata,
                                           target_sum=target_sum,
                                           layer=layer,
                                           copy=True)
        sc.pp.log1p(adata_norm,
                    layer=layer)
        adata_counts = adata_norm.X
        if layer is not None:
            adata_counts = adata_norm.layers[layer]
    if normalize_total and not log1p:
        adata_norm = sc.pp.normalize_total(adata,
                                           target_sum=target_sum,
                                           layer=layer,
                                           copy=True)
        adata_counts = adata_norm.X
        if layer is not None:
            adata_counts = adata_norm.layers[layer]
    if not normalize_total and log1p:
        adata_log1p = sc.pp.log1p(adata,
                                  layer=layer,
                                  copy=True)
        adata_counts = adata_log1p.X
        if layer is not None:
            adata_counts = adata_log1p.layers[layer]
    return adata_counts


def _get_psd(adata,
             gene_id,
             gene_age,
             keep='min',
             layer=None,
             normalize_total=False,
             log1p=False,
             target_sum=1e6):
    """
    A helper function to pre-process AnnData.

    :param adata: AnnData object of shape n_obs × n_vars. Rows correspond to cells and columns to genes.
    :param gene_id: Expects GeneID column from orthomap DataFrame.
    :param gene_age: Expects Phylostratum column from orthomap DataFrame.
    :param keep: In case of duplicated GeneIDs with different Phylostrata assignments, either keep 'min' or 'max' value.
    :param layer: Layer to work on instead of X. If None, X is used.
    :param normalize_total: Normalize counts per cell prior TEI calculation.
    :param log1p: Logarithmize the data matrix prior TEI calculation.
    :param target_sum: After normalization, each observation (cell) has a total count equal to target_sum.
    :return: A list of results such as:
             var_names_df, id_age_df_keep_subset, adata_counts, var_names_subset, sumx, sumx_recd, ps, psd

    :type adata: AnnData
    :type gene_id: list
    :type gene_age: list
    :type keep: str
    :type layer: str
    :type normalize_total: bool
    :type log1p: bool
    :type target_sum: float
    :rtype: list

    Example
    -------
    >>> import scanpy as sc
    >>> from orthomap import orthomap2tei, datasets
    >>> # download pre-calculated orthomap
    >>> #query_orthomap = orthomap2tei.read_orthomap(orthomapfile='Sun2021_Orthomap.tsv')
    >>> sun21_orthomap_file = datasets.sun21_orthomap(datapath='.')
    >>> # load query species orthomap
    >>> query_orthomap = orthomap2tei.read_orthomap(orthomapfile=sun21_orthomap_file)
    >>> # download and load scRNA data
    >>> #packer19_small = sc.read('packer19_small.h5ad')
    >>> packer19_small = datasets.packer19_small(datapath='.')
    >>> # get psd from existing adata object
    >>> celegans_var_names_df,
    >>>     celegans_id_age_df_keep_subset,
    >>>     celegans_adata_counts,
    >>>     celegans_var_names_subset,
    >>>     celegans_sumx,
    >>>     celegans_sumx_recd,
    >>>     celegans_ps,
    >>>     celegans_psd = orthomap2tei._get_psd(
    >>>     adata=packer19_small,
    >>>     gene_id=query_orthomap['GeneID'],
    >>>     gene_age=query_orthomap['Phylostratum'])
    """
    id_age_df = pd.DataFrame(data={'GeneID': gene_id,
                                   'Phylostrata': gene_age})
    # check and drop duplicated GeneID
    id_age_df_keep = _keep_min_max(df=id_age_df,
                                   keep=keep,
                                   dup_col='GeneID',
                                   sort_col='Phylostrata')
    # get overlap with var_names before NaN removal
    var_names_df = pd.merge(left=pd.DataFrame(adata.var_names,
                                              columns=['GeneID']),
                            right=id_age_df_keep,
                            how='left',
                            on='GeneID')
    # check and remove NaN
    id_age_df_keep = id_age_df_keep[~id_age_df_keep['Phylostrata'].isna()]
    # get overlap
    gene_intersection = pd.Index(id_age_df_keep['GeneID']).intersection(adata.var_names)
    id_age_df_keep_subset = id_age_df_keep.loc[id_age_df_keep['GeneID'].isin(gene_intersection)]
    id_age_df_keep_subset = id_age_df_keep_subset.sort_values('GeneID')
    adata_counts = _get_counts(adata=adata,
                               layer=layer,
                               normalize_total=normalize_total,
                               log1p=log1p,
                               target_sum=target_sum)
    adata_counts = adata_counts[:, adata.var_names.isin(id_age_df_keep_subset['GeneID'])]
    var_names_subset = adata.var_names[adata.var_names.isin(id_age_df_keep_subset['GeneID'])]
    var_names_subset, \
        var_names_subset_idx = var_names_subset.sort_values(return_indexer=True)
    adata_counts = adata_counts[:, var_names_subset_idx]
    sumx = adata_counts.sum(1)
    sumx_rec = np.reciprocal(sumx)
    sumx_recd = scipy.sparse.diags(np.array(sumx_rec).flatten())
    ps = np.array(id_age_df_keep_subset['Phylostrata'])
    psd = scipy.sparse.diags(ps)
    return [var_names_df,
            id_age_df_keep_subset,
            adata_counts,
            var_names_subset,
            sumx,
            sumx_recd,
            ps,
            psd]


def add_gene_age2adata_var(adata,
                           gene_id,
                           gene_age,
                           keep='min',
                           var_name='Phylostrata'):
    """
    This function add gene age to an existing AnnData object.

    :param adata: AnnData object of shape n_obs × n_vars. Rows correspond to cells and columns to genes.
    :param gene_id: Expects GeneID column from orthomap DataFrame.
    :param gene_age: Expects GeneID column from orthomap DataFrame.
    :param keep: In case of duplicated GeneIDs with different Phylostrata assignments, either keep 'min' or 'max' value.
    :param var_name: Variable name to be used for gene age values in existing AnnData object.
    :return: Altered AnnData.

    :type adata: AnnData
    :type gene_id: list
    :type gene_age: list
    :type keep: str
    :type var_name: str
    :rtype: AnnData

    Example
    -------
    >>> import scanpy as sc
    >>> from orthomap import datasets, orthomap2tei
    >>> # download pre-calculated orthomap
    >>> #query_orthomap = orthomap2tei.read_orthomap(orthomapfile='Sun2021_Orthomap.tsv')
    >>> sun21_orthomap_file = datasets.sun21_orthomap(datapath='.')
    >>> # load query species orthomap
    >>> query_orthomap = orthomap2tei.read_orthomap(orthomapfile=sun21_orthomap_file)
    >>> # download and load scRNA data
    >>> #packer19_small = sc.read('packer19_small.h5ad')
    >>> packer19_small = datasets.packer19_small(datapath='.')
    >>> # add gene age values to existing adata object
    >>> orthomap2tei.add_gene_age2adata_var(
    >>>     adata=packer19_small,
    >>>     gene_id=query_orthomap['GeneID'],
    >>>     gene_age=query_orthomap['Phylostratum'])
    >>> packer19_small.var
    """
    id_age_df = pd.DataFrame(data={'GeneID': gene_id,
                                   'Phylostrata': gene_age})
    # check and drop duplicated GeneID
    id_age_df_keep = _keep_min_max(df=id_age_df,
                                   keep=keep,
                                   dup_col='GeneID',
                                   sort_col='Phylostrata')
    # get overlap with var_names
    var_names_df = pd.merge(left=pd.DataFrame(adata.var_names.values,
                                              columns=['GeneID']),
                            right=id_age_df_keep,
                            how='left',
                            on='GeneID')
    # add gene age
    adata.var[var_name] = list(var_names_df['Phylostrata'])


def _get_pairwise_comb_self(list1, exclude_self=True):
    """
    A helper function to get all pairwise combinations as a list of tuples.

    :param list1: List to use for combinations.
    :param exclude_self:
    :return: A list of combinations as tuples.

    :type list1: list
    :type exclude_self: bool
    :rtype: list
    """
    pairwise_comb = []
    for l1_idx, l1 in enumerate(list1[:-1]):
        for l2 in list1[l1_idx:]:
            if exclude_self:
                if l1 != l2:
                    pairwise_comb.append((l1, l2))
            else:
                pairwise_comb.append((l1, l2))
    if not exclude_self:
        pairwise_comb.append((list1[-1], list1[-1]))
    return pairwise_comb


def _get_pairwise_comb(list1, list2):
    """
    A helper function to get all pairwise combinations from two lists.

    :param list1: List1 to use for combinations.
    :param list2: List1 to use for combinations.
    :return: A list of combinations as tuples.

    :type list1: list
    :type list2: list
    :rtype: list
    """
    pairwise_comb = []
    for l1 in list1:
        for l2 in list2:
            pairwise_comb.append((l1, l2))
    return pairwise_comb


def get_tei(adata,
            gene_id,
            gene_age,
            keep='min',
            layer=None,
            add_var=True,
            var_name='Phylostrata',
            add_obs=True,
            obs_name='tei',
            boot=False,
            bt=10,
            normalize_total=False,
            log1p=False,
            target_sum=1e6):
    """
    This function computes the phylogenetically based transcriptome evolutionary
    index (TEI) similar to Domazet-Loso & Tautz, 2010.

    The TEI measure represents the weighted arithmetic mean
    (expression levels as weights for the phylostratum value) over all
    evolutionary age categories denoted as _phylostra_.

    :: math::
        \deqn{TEI_s = sum (e_is * ps_i) / sum e_is}

    where \eqn{TEI_s} denotes the TEI value in developmental stage \eqn{s, e_is}
    denotes the gene expression level of gene \eqn{i} in stage \eqn{s}, and \eqn{ps_i}
    denotes the corresponding phylostratum of gene \eqn{i, i = 1,...,N} and
    \eqn{N = total number of genes}.

    If the parameter boot is set to true,
    the strata values are sampled and the global TEI
    is calculated bt times.

    :param adata: AnnData object of shape n_obs × n_vars. Rows correspond to cells and columns to genes.
    :param gene_id: Expects GeneID column from orthomap DataFrame.
    :param gene_age: Expects GeneID column from orthomap DataFrame.
    :param keep: In case of duplicated GeneIDs with different Phylostrata assignments, either keep 'min' or 'max' value.
    :param layer: Layer to work on instead of X. If None, X is used.
    :param add_var: Add gene age values as variable to existing AnnData object using var_name.
    :param var_name: Variable name to be used for gene age values in existing AnnData object.
    :param add_obs: Add TEI values as observation to existing AnnData object using obs_name.
    :param obs_name: Observation name to be used for TEI values in existing AnnData object.
    :param boot: Specify if bootstrap TEI values should be calculated and returned as DataFrame.
    :param bt: Number of bootstrap to calculate.
    :param normalize_total: Normalize counts per cell prior TEI calculation.
    :param log1p: Logarithmize the data matrix prior TEI calculation.
    :param target_sum: After normalization, each observation (cell) has a total count equal to target_sum.
    :return: Transcriptome evolutionary index (TEI) values.

    :type adata: AnnData
    :type gene_id: list
    :type gene_age: list
    :type keep: str
    :type layer: str
    :type add_var: bool
    :type var_name: str
    :type add_obs: bool
    :type obs_name: str
    :type boot: bool
    :type bt: int
    :type normalize_total: bool
    :type log1p: bool
    :type target_sum: float
    :rtype: pandas.DataFrame

    Example
    -------
    >>> import scanpy as sc
    >>> import matplotlib.pyplot as plt
    >>> import seaborn as sns
    >>> from statannot import add_stat_annotation
    >>> from orthomap import datasets, orthomap2tei
    >>> # download pre-calculated orthomap
    >>> #query_orthomap = orthomap2tei.read_orthomap(orthomapfile='Sun2021_Orthomap.tsv')
    >>> sun21_orthomap_file = datasets.sun21_orthomap(datapath='.')
    >>> # load query species orthomap
    >>> query_orthomap = orthomap2tei.read_orthomap(orthomapfile=sun21_orthomap_file)
    >>> # download and load scRNA data
    >>> #packer19_small = sc.read('packer19_small.h5ad')
    >>> packer19_small = datasets.packer19_small(datapath='.')
    >>> # add TEI values to existing adata object
    >>> orthomap2tei.get_tei(
    >>>     adata=packer19_small,
    >>>     gene_id=query_orthomap['GeneID'],
    >>>     gene_age=query_orthomap['Phylostratum'],
    >>>     add_var=True,
    >>>     add_obs=True)
    >>> # plot tei boxplot grouped by embryo.time.bin observation
    >>> ax = sns.boxplot(
    >>>     x='embryo.time.bin',
    >>>     y='tei',
    >>>     data=packer19_small.obs)
    >>> test_results = add_stat_annotation(
    >>>     ax,
    >>>     x='embryo.time.bin',
    >>>     y='tei',
    >>>     data=packer19_small.obs,
    >>>     box_pairs=orthomap2tei._get_pairwise_comb_self(
    >>>         list1=packer19_small.obs['embryo.time.bin'].value_counts().index),
    >>>     test='Mann-Whitney',
    >>>     text_format='star',
    >>>     loc='outside',
    >>>     verbose=2)
    >>> plt.show()
    >>> # plot tei violinplot for each cell.type grouped by cell.type and embryo.time.bin observation
    >>> # create new observation as a combination from embryo.time.bin and cell.type
    >>> packer19_small.obs['etb_cell.type'] = packer19_small.obs[['embryo.time.bin', 'cell.type']].apply(
    >>>     lambda x: str(x[0]) + '_' + x[1], axis=1)
    >>> # convert into category
    >>> packer19_small.obs['etb_cell.type'] = packer19_small.obs['etb_cell.type'].astype('category')
    >>> # reorder categories
    >>> packer19_small.obs['etb_cell.type'] = packer19_small.obs['etb_cell.type'].cat
    >>>     .reorder_categories(list(packer19_small.obs['etb_cell.type']
    >>>     .value_counts().index[np.argsort([int(x.split('_')[0]) for x in
    >>>     list(packer19_small.obs['etb_cell.type'].value_counts().index)])]))
    >>> for c in packer19_small.obs['cell.type'].value_counts().index:
    >>>     plt.figure()
    >>>     sns.violinplot(
    >>>     x=packer19_small.obs[packer19_small.obs['cell.type'].isin([c])]['etb_cell.type'].cat
    >>>     .remove_unused_categories(),
    >>>     y='tei',
    >>>     data=packer19_small.obs[packer19_small.obs['cell.type'].isin([c])])
    >>> plt.show()
    >>> # get 10 bootstrap TEI values
    >>> orthomap2tei.get_tei(
    >>>     adata=packer19_small,
    >>>     gene_id=query_orthomap['GeneID'],
    >>>     gene_age=query_orthomap['Phylostratum'],
    >>>     boot=True,
    >>>     bt=10)
    """
    var_names_df,\
        id_age_df_keep_subset,\
        adata_counts,\
        var_names_subset,\
        sumx,\
        sumx_recd,\
        ps,\
        psd = _get_psd(adata=adata,
                       gene_id=gene_id,
                       gene_age=gene_age,
                       keep=keep,
                       layer=layer,
                       normalize_total=normalize_total,
                       log1p=log1p,
                       target_sum=target_sum)
    wmatrix = psd.dot(adata_counts.transpose()).transpose()
    pmatrix = sumx_recd.dot(wmatrix)
    tei = pmatrix.sum(1)
    tei_df = pd.DataFrame(tei,
                          columns=[obs_name])
    tei_df.index = adata.obs_names
    if add_var:
        add_gene_age2adata_var(adata=adata,
                               gene_id=gene_id,
                               gene_age=gene_age,
                               keep=keep,
                               var_name=var_name)
    if add_obs:
        adata.obs[obs_name] = tei_df
    if boot:
        with alive_bar(bt) as bar:
            for i in range(bt):
                np.random.shuffle(ps)
                psd = scipy.sparse.diags(ps)
                tei = sumx_recd.dot(psd.dot(adata_counts.transpose()).transpose()).sum(1)
                tei_df[i] = tei
                bar()
    return tei_df


def get_pmatrix(adata,
                gene_id,
                gene_age,
                keep='min',
                layer=None,
                layer_name='pmatrix',
                add_var=True,
                add_obs=True,
                normalize_total=False,
                log1p=False,
                target_sum=1e6):
    """
    This function computes the partial transcriptome evolutionary index (TEI) values for each single gene.

    Prior TEI calculation, counts can be normalized (default: False) to a total count number (default: 1e6) and
    log transformed (default: False).

    In detail, each gene gets a TEI contribution profile as follows:
    \deqn{TEI_is = f_is * ps_i}
    where TEI_is is the partial TEI value of gene i,
    \eqn{f_is = e_is / \sum e_is} and \eqn{ps_i} is the phylostratum of gene i.

    The partial TEI matrix can be used to perform different cluster
    analyses and also gives an overall impression of the contribution of each
    gene to the global TEI pattern.

    :param adata: AnnData object of shape n_obs × n_vars. Rows correspond to cells and columns to genes.
    :param gene_id: Expects GeneID column from orthomap DataFrame.
    :param gene_age: Expects Phylostratum column from orthomap DataFrame.
    :param keep: Either define 'min' (ascending pre-sorting) or 'max' (non-ascending pre-sorting) to keep duplicates.
    :param layer: Layer to work on instead of X. If None, X is used.
    :param layer_name: Layer to add to existing AnnData object.
    :param add_var: Add original variables to new AnnData object.
    :param add_obs: Add original observations to new AnnData object.
    :param normalize_total: Normalize counts per cell prior TEI calculation.
    :param log1p: Logarithmize the data matrix prior TEI calculation.
    :param target_sum: After normalization, each observation (cell) has a total count equal to target_sum.
    :return: Partial transcriptome evolutionary index (TEI) values.

    :type adata: AnnData
    :type gene_id: list
    :type gene_age: list
    :type keep: str
    :type layer: str
    :type layer_name: str
    :type add_var: bool
    :type add_obs: bool
    :type normalize_total: bool
    :type log1p: bool
    :type target_sum: float
    :rtype: AnnData

    Example
    -------
    >>> import scanpy as sc
    >>> from orthomap import orthomap2tei, datasets
    >>> # download pre-calculated orthomap
    >>> #query_orthomap = orthomap2tei.read_orthomap(orthomapfile='Sun2021_Orthomap.tsv')
    >>> sun21_orthomap_file = datasets.sun21_orthomap(datapath='.')
    >>> # load query species orthomap
    >>> query_orthomap = orthomap2tei.read_orthomap(orthomapfile=sun21_orthomap_file)
    >>> # download and load scRNA data
    >>> #packer19_small = sc.read('packer19_small.h5ad')
    >>> packer19_small = datasets.packer19_small(datapath='.')
    >>> # get pmatrix as new adata object
    >>> packer19_small_pmatrix = orthomap2tei.get_pmatrix(
    >>>     adata=packer19_small,
    >>>     gene_id=query_orthomap['GeneID'],
    >>>     gene_age=query_orthomap['Phylostratum'])
    """
    var_names_df,\
        id_age_df_keep_subset,\
        adata_counts,\
        var_names_subset,\
        sumx,\
        sumx_recd,\
        ps,\
        psd = _get_psd(adata=adata,
                       gene_id=gene_id,
                       gene_age=gene_age,
                       keep=keep,
                       layer=layer,
                       normalize_total=normalize_total,
                       log1p=log1p,
                       target_sum=target_sum)
    wmatrix = psd.dot(adata_counts.transpose()).transpose()
    pmatrix = sumx_recd.dot(wmatrix)
    adata_pmatrix = ad.AnnData(adata_counts)
    adata_pmatrix.layers[layer_name] = pmatrix
    adata_pmatrix.obs_names = adata.obs_names
    adata_pmatrix.var_names = var_names_subset
    if add_obs:
        for ko in adata.obs.keys():
            adata_pmatrix.obs[ko] = adata.obs[ko]
    if add_var:
        for kv in adata.var.keys():
            adata_pmatrix.var[kv] = pd.merge(left=adata_pmatrix.var,
                                             right=adata.var[kv][adata.var_names.isin(id_age_df_keep_subset['GeneID'])],
                                             left_index=True,
                                             right_index=True)[kv]
    adata_pmatrix.var['Phylostrata'] = list(pd.merge(left=pd.DataFrame(adata_pmatrix.var_names,
                                                                       columns=['GeneID']),
                                                     right=var_names_df,
                                                     how='left',
                                                     on='GeneID')['Phylostrata'])
    return adata_pmatrix


def get_pstrata(adata,
                gene_id,
                gene_age,
                keep='min',
                layer=None,
                cumsum=False,
                group_by_obs=None,
                obs_fillna='__NaN',
                obs_type='mean',
                standard_scale=None,
                normalize_total=False,
                log1p=False,
                target_sum=1e6):
    """
    This function computes the partial transcriptome evolutionary index (TEI) values combined for each stratum.

    The resulting values can be combined per observation group e.g.: pre-defined cell types (default: None),
    according to the selected observation type (default:'mean') and further scaled between 0 and 1 (default: None)
    either per var (standard_scale=0) or per obs (standard_scale=1).

    Prior TEI calculation, counts can be normalized (default: False) to a total count number (default: 1e6) and
    log transformed (default: False).

    In detail, each gene gets a TEI contribution profile as follows:
    \deqn{TEI_is = f_is * ps_i}
    where TEI_is is the partial TEI value of gene i,
    \eqn{f_is = e_is / \sum e_is} and \eqn{ps_i} is the phylostratum of gene i.

    \eqn{TEI_is} values are combined per \eqn{ps}.

    The partial TEI values combined per strata give an overall impression
    of the contribution of each stratum to the global TEI pattern.

    :param adata: AnnData object of shape n_obs × n_vars. Rows correspond to cells and columns to genes.
    :param gene_id: Expects GeneID column from orthomap DataFrame.
    :param gene_age: Expects Phylostratum column from orthomap DataFrame.
    :param keep: Either define 'min' (ascending pre-sorting) or 'max' (non-ascending pre-sorting) to keep duplicates.
    :param layer: Layer to work on instead of X. If None, X is used.
    :param cumsum: Return cumsum.
    :param group_by_obs: AnnData observation to be used as a group to combine partial transcriptome evolutionary index
                         (TEI) values.
    :param obs_fillna: Specify how NaN values should be named for observation.
    :param obs_type: Specify how values should be combined per observation group. Possible values are 'mean', 'median',
                     'sum', 'min' and 'max'.
    :param standard_scale: Wether or not to standardize the given axis (0: colums, 1: rows) between 0 and 1,
                           meaning for each variable or group, subtract the minimum and divide each by its maximum.
    :param normalize_total: Normalize counts per cell prior TEI calculation.
    :param log1p: Logarithmize the data matrix prior TEI calculation.
    :param target_sum: After normalization, each observation (cell) has a total count equal to target_sum.
    :return: List of two DataFrame. First DataFrame contains the summed partial TEI values per strata.
             Second DataFrame contains summed partial TEI values divided by the corresponding global TEI value,
             which represent percentage of global TEI per strata.

    :type adata: AnnData
    :type gene_id: list
    :type gene_age: list
    :type keep: str
    :type layer: str
    :type cumsum: bool
    :type group_by_obs: str
    :type obs_fillna: str
    :type obs_type: str
    :type standard_scale: int
    :type normalize_total: bool
    :type log1p: bool
    :type target_sum: float
    :rtype: list

    Example
    -------
    >>> import scanpy as sc
    >>> import matplotlib.pyplot as plt
    >>> import seaborn as sns
    >>> from orthomap import orthomap2tei, datasets
    >>> # download pre-calculated orthomap
    >>> #query_orthomap = orthomap2tei.read_orthomap(orthomapfile='Sun2021_Orthomap.tsv')
    >>> sun21_orthomap_file = datasets.sun21_orthomap(datapath='.')
    >>> # load query species orthomap
    >>> query_orthomap = orthomap2tei.read_orthomap(orthomapfile=sun21_orthomap_file)
    >>> # download and load scRNA data
    >>> #packer19_small = sc.read('packer19_small.h5ad')
    >>> packer19_small = datasets.packer19_small(datapath='.')
    >>> # get pstrata
    >>> packer19_small_pstrata = orthomap2tei.get_pstrata(
    >>>     adata=packer19_small,
    >>>     gene_id=query_orthomap['GeneID'],
    >>>     gene_age=query_orthomap['Phylostratum'])
    >>> # get cumsum over strata
    >>> packer19_small_pstrata_cumsum = orthomap2tei.get_pstrata(
    >>>     adata=packer19_small,
    >>>     gene_id=query_orthomap['GeneID'],
    >>>     gene_age=query_orthomap['Phylostratum'],
    >>>     cumsum=True)
    >>> # group by embryo.time.bin observation
    >>> packer19_small_pstrata_grouped = orthomap2tei.get_pstrata(
    >>>     adata=packer19_small,
    >>>     gene_id=query_orthomap['GeneID'],
    >>>     gene_age=query_orthomap['Phylostratum'],
    >>>     group_by_obs='embryo.time.bin')
    >>> # plot strata as lines
    >>> packer19_small_pstrata_grouped[0].transpose().plot.line(stacked=True, cmap='Accent')
    >>> # plot heatmap using partial TEI values
    >>> sns.heatmap(packer19_small_pstrata_grouped[0], annot=True, cmap='viridis')
    >>> plt.show()
    >>> # plot heatmap using partial TEI percent
    >>> sns.heatmap(packer19_small_pstrata_grouped[1], annot=True, cmap='viridis')
    >>> plt.show()
    """
    var_names_df,\
        id_age_df_keep_subset,\
        adata_counts,\
        var_names_subset,\
        sumx,\
        sumx_recd,\
        ps,\
        psd = _get_psd(adata=adata,
                       gene_id=gene_id,
                       gene_age=gene_age,
                       keep=keep,
                       layer=layer,
                       normalize_total=normalize_total,
                       log1p=log1p,
                       target_sum=target_sum)
    wmatrix = psd.dot(adata_counts.transpose()).transpose()
    pmatrix = sumx_recd.dot(wmatrix)
    tei = pmatrix.sum(1)
    phylostrata = list(set(id_age_df_keep_subset['Phylostrata']))
    pstrata_norm_by_sumx = np.zeros((len(phylostrata), pmatrix.shape[0]))
    pstrata_norm_by_pmatrix_sum = np.zeros((len(phylostrata), pmatrix.shape[0]))
    for pk_idx, pk in enumerate(phylostrata):
        pstrata_norm_by_sumx[pk_idx, ] = np.array(pmatrix[:, id_age_df_keep_subset['Phylostrata'].isin([pk])]
                                                  .sum(1)).flatten()
        pstrata_norm_by_pmatrix_sum[pk_idx, ] = np.array(pmatrix[:, id_age_df_keep_subset['Phylostrata']
                                                         .isin([pk])].sum(1) / tei).flatten()
    pstrata_norm_by_sumx_df = pd.DataFrame(pstrata_norm_by_sumx)
    pstrata_norm_by_sumx_df['ps'] = phylostrata
    pstrata_norm_by_sumx_df.set_index('ps',
                                      inplace=True)
    pstrata_norm_by_sumx_df.columns = adata.obs_names
    pstrata_norm_by_pmatrix_sum_df = pd.DataFrame(pstrata_norm_by_pmatrix_sum)
    pstrata_norm_by_pmatrix_sum_df['ps'] = phylostrata
    pstrata_norm_by_pmatrix_sum_df.set_index('ps',
                                             inplace=True)
    pstrata_norm_by_pmatrix_sum_df.columns = adata.obs_names
    if cumsum:
        pstrata_norm_by_sumx_df = pstrata_norm_by_sumx_df.cumsum(0)
        pstrata_norm_by_pmatrix_sum_df = pstrata_norm_by_pmatrix_sum_df.cumsum(0)
    if group_by_obs is not None:
        if adata.obs[group_by_obs].dtype.name == 'category':
            obs_group_nan = pd.DataFrame(adata.obs[group_by_obs])
        else:
            obs_group_nan = pd.DataFrame(adata.obs[group_by_obs].fillna(obs_fillna))
        if obs_type == 'mean':
            pstrata_norm_by_sumx_df =\
                pstrata_norm_by_sumx_df.transpose().groupby(obs_group_nan[group_by_obs]).mean().transpose()
            pstrata_norm_by_pmatrix_sum_df =\
                pstrata_norm_by_pmatrix_sum_df.transpose().groupby(obs_group_nan[group_by_obs]).mean().transpose()
        if obs_type == 'median':
            pstrata_norm_by_sumx_df =\
                pstrata_norm_by_sumx_df.transpose().groupby(obs_group_nan[group_by_obs]).median().transpose()
            pstrata_norm_by_pmatrix_sum_df =\
                pstrata_norm_by_pmatrix_sum_df.transpose().groupby(obs_group_nan[group_by_obs]).median().transpose()
        if obs_type == 'sum':
            pstrata_norm_by_sumx_df =\
                pstrata_norm_by_sumx_df.transpose().groupby(obs_group_nan[group_by_obs]).sum().transpose()
            pstrata_norm_by_pmatrix_sum_df =\
                pstrata_norm_by_pmatrix_sum_df.transpose().groupby(obs_group_nan[group_by_obs]).sum().transpose()
        if obs_type == 'min':
            pstrata_norm_by_sumx_df =\
                pstrata_norm_by_sumx_df.transpose().groupby(obs_group_nan[group_by_obs]).min().transpose()
            pstrata_norm_by_pmatrix_sum_df =\
                pstrata_norm_by_pmatrix_sum_df.transpose().groupby(obs_group_nan[group_by_obs]).min().transpose()
        if obs_type == 'max':
            pstrata_norm_by_sumx_df =\
                pstrata_norm_by_sumx_df.transpose().groupby(obs_group_nan[group_by_obs]).max().transpose()
            pstrata_norm_by_pmatrix_sum_df =\
                pstrata_norm_by_pmatrix_sum_df.transpose().groupby(obs_group_nan[group_by_obs]).max().transpose()
    if standard_scale is not None:
        if standard_scale == 0:
            pstrata_norm_by_sumx_df = pstrata_norm_by_sumx_df.apply(_min_max_to_01,
                                                                    axis=1,
                                                                    raw=True)
            pstrata_norm_by_pmatrix_sum_df = pstrata_norm_by_pmatrix_sum_df.apply(_min_max_to_01,
                                                                                  axis=1,
                                                                                  raw=True)
        if standard_scale == 1:
            pstrata_norm_by_sumx_df = pstrata_norm_by_sumx_df.apply(_min_max_to_01,
                                                                    axis=0,
                                                                    raw=True)
            pstrata_norm_by_pmatrix_sum_df = pstrata_norm_by_pmatrix_sum_df.apply(_min_max_to_01,
                                                                                  axis=0,
                                                                                  raw=True)
    return [pstrata_norm_by_sumx_df,
            pstrata_norm_by_pmatrix_sum_df]


def _min_max_to_01(ndarray):
    """
    A helper function to standardize data between 0 and 1,
    meaning for each variable, subtract the minimum and divide each by its maximum.

    :param ndarray: Data to be standardized.
    :return: Standardized data.

    :type ndarray: numpy.ndarray
    :rtype: numpy.ndarray

    Example
    -------
    >>> import numpy as np
    >>> from orthomap import orthomap2tei
    >>> random_array = np.random.rand(10)
    >>> orthomap2tei._min_max_to_01(random_array)
    """
    ndarray_min = np.min(ndarray)
    ndarray_max = np.max(ndarray)
    if ndarray_min == ndarray_max:
        ndarray_min_max = [x - ndarray_min for x in ndarray]
    else:
        ndarray_min_max = [((x - ndarray_min)/(ndarray_max-ndarray_min)) for x in ndarray]
    return ndarray_min_max


def get_ematrix(adata,
                layer=None,
                group_by_var=None,
                var_type='mean',
                var_fillna='__NaN',
                group_by_obs=None,
                obs_type='mean',
                obs_fillna='__NaN',
                standard_scale=None,
                normalize_total=False,
                log1p=False,
                target_sum=1e6):
    """
    This function computes expression profiles for all genes or group of genes 'group_by_var' (default: None).

    The expression values are first combined per var type 'var_type' (default: mean).

    The resulting values can be combined per observation group 'group_by_obs' e.g.: pre-defined cell types
    (default: None), according to the selected observation type 'obs_type' (default:'mean') and further scaled between
    0 and 1 (default: None) either per var (standard_scale=0) or per obs (standard_scale=1).

    In detail, if standard_scale axis is set to None, the var_type mean/median/sum expression is being computed over
    cells and, if group_by_obs is not None, combined per given obs group by mean/median/sum.

    In detail, if standard_scale axis is set to 0, the mean/median/sum relative expression profile is being computed
    over cells and, if group_by_obs is not None, combined per given obs group by mean/median/sum as follows:

    f_c = (e_c - e_min)/(e_max - e_min)

    where e_min and e_max denote either the minimum/maximum mean/median/sum
    expression level over cells c.

    In detail, if standard_scale axis is set to 1, the mean/median/sum relative expression profile is being computed
    over gene age classes (phylostrata) and, if group_by_obs is not None, combined per given obs group by
    mean/median/sum as follows:

    f_ps = (e_ps - e_min)/(e_max - e_min)

    where e_min and e_max denote either the minimum/maximum mean/median/sum
    expression level over gene age class (phylostrata ps).

    This linear transformation corresponds to a shift by e_min -
    e_max. As a result, the relative expression level f_c of cell c or f_ps
    of phylotstratum ps with minimum e_c or e_ps is 0,
    whereas the relative expression level f_c of cell c or f_ps of phylotstratum ps
    with maximum e_c or e_ps is 1, and the relative
    expression levels of all other cells c or
    phylostrata ps range between 0 and 1.

    :param adata: AnnData object of shape n_obs × n_vars. Rows correspond to cells and columns to genes.
    :param layer: Layer to work on instead of X. If None, X is used.
    :param group_by_var: AnnData variable to be used as a group to combine count values.
    :param var_type: Specify how values should be combined per variable group. Possible values are 'mean', 'median',
                     'sum', 'min' and 'max'.
    :param var_fillna: Specify how NaN values should be named for variable.
    :param group_by_obs: AnnData observation to be used as a group to combine count values.
    :param obs_type: Specify how values should be combined per observation group. Possible values are 'mean', 'median',
                     'sum', 'min' and 'max'.
    :param obs_fillna: Specify how NaN values should be named for observation.
    :param standard_scale: Wether or not to standardize the given axis (0: colums, 1: rows) between 0 and 1,
                           meaning for each variable or group, subtract the minimum and divide each by its maximum.
    :param normalize_total: Normalize counts per cell.
    :param log1p: Logarithmize the data matrix.
    :param target_sum: After normalization, each observation (cell) has a total count equal to target_sum.
    :return: Expression profile DataFrame.

    :type adata: AnnData
    :type layer: str
    :type group_by_var: str
    :type var_type: str
    :type var_fillna: str
    :type group_by_obs: str
    :type obs_type: str
    :type obs_fillna: str
    :type standard_scale: int
    :type normalize_total: bool
    :type log1p: bool
    :type target_sum: float
    :rtype: pandas.DataFrame

    Example
    -------
    >>> import scanpy as sc
    >>> import matplotlib.pyplot as plt
    >>> import seaborn as sns
    >>> from orthomap import orthomap2tei, datasets
    >>> # download pre-calculated orthomap
    >>> #query_orthomap = orthomap2tei.read_orthomap(orthomapfile='Sun2021_Orthomap.tsv')
    >>> sun21_orthomap_file = datasets.sun21_orthomap(datapath='.')
    >>> # load query species orthomap
    >>> query_orthomap = orthomap2tei.read_orthomap(orthomapfile=sun21_orthomap_file)
    >>> # download and load scRNA data
    >>> #packer19_small = sc.read('packer19_small.h5ad')
    >>> packer19_small = datasets.packer19_small(datapath='.')
    >>> # add gene age values to existing adata object
    >>> orthomap2tei.add_gene_age2adata_var(
    >>>     adata=packer19_small,
    >>>     gene_id=query_orthomap['GeneID'],
    >>>     gene_age=query_orthomap['Phylostratum'])
    >>> packer19_small_ematrix_grouped = orthomap2tei.get_ematrix(
    >>>     adata=packer19_small,
    >>>     group_by_var='Phylostrata',
    >>>     group_by_obs='embryo.time.bin')
    >>> packer19_small_ematrix_grouped.transpose().plot.line(stacked=True, cmap='Accent')
    >>> plt.show()
    >>> sns.heatmap(packer19_small_ematrix_grouped, annot=True, cmap='viridis')
    >>> plt.show()
    >>> # normalize counts (transcript per million - tpm)
    >>> packer19_small_ematrix_grouped_tpm = orthomap2tei.get_ematrix(
    >>>     adata=packer19_small,
    >>>     group_by_var='Phylostrata',
    >>>     group_by_obs='embryo.time.bin',
    >>>     normalize_total=True)
    >>> packer19_small_ematrix_grouped_tpm.transpose().plot.line(stacked=True, cmap='Accent')
    >>> plt.show()
    >>> sns.heatmap(packer19_small_ematrix_grouped_tpm, annot=True, cmap='viridis')
    >>> plt.show()
    """
    adata_counts = _get_counts(adata=adata,
                               layer=layer,
                               normalize_total=normalize_total,
                               log1p=log1p,
                               target_sum=target_sum)
    if group_by_var is not None:
        if adata.var[group_by_var].dtype.name == 'category':
            var_grouped = pd.DataFrame(adata.var[group_by_var]).groupby(group_by_var)
        else:
            var_grouped = pd.DataFrame(adata.var[group_by_var].fillna(var_fillna)).groupby(group_by_var)
        var_groups = var_grouped.groups.keys()
    if group_by_obs is not None:
        if adata.obs[group_by_obs].dtype.name == 'category':
            obs_grouped = pd.DataFrame(adata.obs[group_by_obs]).groupby(group_by_obs)
        else:
            obs_grouped = pd.DataFrame(adata.obs[group_by_obs].fillna(obs_fillna)).groupby(group_by_obs)
        obs_groups = obs_grouped.groups.keys()
    if group_by_var is not None and group_by_obs is not None:
        ematrix_df = pd.DataFrame(np.zeros((len(var_groups), len(obs_groups)), dtype=np.float64),
                                  columns=list(obs_groups),
                                  index=list(var_groups)
                                  )
        for var_idx, (var_group, var_group_idx) in enumerate(var_grouped.indices.items()):
            for obs_idx, (obs_group, obs_group_idx) in enumerate(obs_grouped.indices.items()):
                adata_counts_x = adata_counts[:, var_group_idx][obs_group_idx, :]
                if var_type == obs_type:
                    if var_type == 'mean':
                        ematrix_df.loc[var_group, obs_group] = adata_counts_x.mean()
                    if var_type == 'median':
                        ematrix_df.loc[var_group, obs_group] = np.median(np.apply_along_axis(
                            np.median, 1, adata_counts_x.toarray()))
                    if var_type == 'sum':
                        ematrix_df.loc[var_group, obs_group] = adata_counts_x.sum()
                    if var_type == 'min':
                        ematrix_df.loc[var_group, obs_group] = adata_counts_x.min()
                    if var_type == 'max':
                        ematrix_df.loc[var_group, obs_group] = adata_counts_x.max()
                else:
                    if var_type == 'mean':
                        if obs_type == 'median':
                            ematrix_df.loc[var_group, obs_group] = np.median(np.ravel(adata_counts_x.mean(1)))
                        if obs_type == 'sum':
                            ematrix_df.loc[var_group, obs_group] = adata_counts_x.mean(1).sum()
                        if obs_type == 'min':
                            ematrix_df.loc[var_group, obs_group] = adata_counts_x.mean(1).min()
                        if obs_type == 'max':
                            ematrix_df.loc[var_group, obs_group] = adata_counts_x.mean(1).max()
                    if var_type == 'median':
                        if obs_type == 'mean':
                            ematrix_df.loc[var_group, obs_group] = np.apply_along_axis(
                                np.median, 1, adata_counts_x.toarray()).mean()
                        if obs_type == 'sum':
                            ematrix_df.loc[var_group, obs_group] = np.apply_along_axis(
                                np.median, 1, adata_counts_x.toarray()).sum()
                        if obs_type == 'min':
                            ematrix_df.loc[var_group, obs_group] = np.apply_along_axis(
                                np.median, 1, adata_counts_x.toarray()).min()
                        if obs_type == 'max':
                            ematrix_df.loc[var_group, obs_group] = np.apply_along_axis(
                                np.median, 1, adata_counts_x.toarray()).max()
                    if var_type == 'sum':
                        if obs_type == 'mean':
                            ematrix_df.loc[var_group, obs_group] = adata_counts_x.sum(1).mean()
                        if obs_type == 'median':
                            ematrix_df.loc[var_group, obs_group] = np.median(np.ravel(adata_counts_x.sum(1)))
                        if obs_type == 'min':
                            ematrix_df.loc[var_group, obs_group] = adata_counts_x.sum(1).min()
                        if obs_type == 'max':
                            ematrix_df.loc[var_group, obs_group] = adata_counts_x.sum(1).max()
                    if var_type == 'min':
                        if obs_type == 'mean':
                            ematrix_df.loc[var_group, obs_group] = adata_counts_x.min(1).mean()
                        if obs_type == 'median':
                            ematrix_df.loc[var_group, obs_group] = np.median(np.ravel(adata_counts_x.min(1)))
                        if obs_type == 'sum':
                            ematrix_df.loc[var_group, obs_group] = adata_counts_x.min(1).sum()
                        if obs_type == 'max':
                            ematrix_df.loc[var_group, obs_group] = adata_counts_x.min(1).max()
                    if var_type == 'max':
                        if obs_type == 'mean':
                            ematrix_df.loc[var_group, obs_group] = adata_counts_x.max(1).mean()
                        if obs_type == 'median':
                            ematrix_df.loc[var_group, obs_group] = np.median(np.ravel(adata_counts_x.max(1)))
                        if obs_type == 'sum':
                            ematrix_df.loc[var_group, obs_group] = adata_counts_x.max(1).sum()
                        if obs_type == 'min':
                            ematrix_df.loc[var_group, obs_group] = adata_counts_x.max(1).min()
    if group_by_var is None and group_by_obs is not None:
        ematrix_df = pd.DataFrame(np.zeros((adata.shape[1], len(obs_groups)), dtype=np.float64),
                                  columns=list(obs_groups),
                                  index=adata.var.index
                                  )
        for obs_idx, (obs_group, obs_group_idx) in enumerate(obs_grouped.indices.items()):
            adata_counts_x = adata_counts[obs_group_idx, :]
            if obs_type == 'mean':
                ematrix_df.loc[:, obs_group] = np.ravel(adata_counts_x.mean(0))
            if obs_type == 'median':
                ematrix_df.loc[:, obs_group] = np.apply_along_axis(np.median, 0, adata_counts_x.toarray())
            if obs_type == 'sum':
                ematrix_df.loc[:, obs_group] = np.ravel(adata_counts_x.sum(0))
            if obs_type == 'min':
                ematrix_df.loc[:, obs_group] = np.ravel(adata_counts_x.min(0))
            if obs_type == 'max':
                ematrix_df.loc[:, obs_group] = np.ravel(adata_counts_x.max(0))
    if group_by_var is not None and group_by_obs is None:
        ematrix_df = pd.DataFrame(np.zeros((len(var_groups), adata.shape[0]), dtype=np.float64),
                                  columns=adata.obs.index,
                                  index=list(var_groups)
                                  )
        for var_idx, (var_group, var_group_idx) in enumerate(var_grouped.indices.items()):
            adata_counts_x = adata_counts[:, var_group_idx]
            if var_type == 'mean':
                ematrix_df.loc[var_group, :] = np.ravel(adata_counts_x.mean(1))
            if var_type == 'median':
                ematrix_df.loc[var_group, :] = np.apply_along_axis(np.median, 1, adata_counts_x.toarray())
            if var_type == 'sum':
                ematrix_df.loc[var_group, :] = np.ravel(adata_counts_x.sum(1))
            if var_type == 'min':
                ematrix_df.loc[var_group, :] = np.ravel(adata_counts_x.min(1))
            if var_type == 'max':
                ematrix_df.loc[var_group, :] = np.ravel(adata_counts_x.max(1))
    if group_by_var is None and group_by_obs is None:
        ematrix_df = pd.DataFrame(adata_counts.toarray().transpose(),
                                  columns=adata.obs.index,
                                  index=adata.var.index
                                  )
    if standard_scale is not None:
        if standard_scale == 0:
            ematrix_df = ematrix_df.apply(_min_max_to_01,
                                          axis=1,
                                          raw=True)
        if standard_scale == 1:
            ematrix_df = ematrix_df.apply(_min_max_to_01,
                                          axis=0,
                                          raw=True)
    return ematrix_df


def get_rematrix(adata,
                 gene_id,
                 gene_age,
                 keep='min',
                 layer=None,
                 use='counts',
                 var_type='mean',
                 group_by_obs=None,
                 obs_fillna='__NaN',
                 obs_type='mean',
                 standard_scale=None,
                 normalize_total=False,
                 log1p=False,
                 target_sum=1e6):
    """
    This function computes relative expression profiles.

    In detail, if standard_scale axis is set to None, the mean/median/sum expression is being computed over cells
    and, if group_by_obs is not None, combined per given obs group by mean/median/sum.

    In detail, if standard_scale axis is set to 0, the mean/median/sum relative expression profile is being computed
    over cells and, if group_by_obs is not None, combined per given obs group by mean/median/sum as follows:

    f_c = (e_c - e_min)/(e_max - e_min)

    where e_min and e_max denote either the minimum/maximum mean/median/sum
    expression level over cells c.

    In detail, if standard_scale axis is set to 1, the mean/median/sum relative expression profile is being computed
    over gene age classes (phylostrata) and, if group_by_obs is not None, combined per given obs group by
    mean/median/sum as follows:

    f_ps = (e_ps - e_min)/(e_max - e_min)

    where e_min and e_max denote either the minimum/maximum mean/median/sum
    expression level over gene age class (phylostrata ps).

    This linear transformation corresponds to a shift by e_min -
    e_max. As a result, the relative expression level f_c of cell c or f_ps
    of phylotstratum ps with minimum e_c or e_ps is 0,
    whereas the relative expression level f_c of cell c or f_ps of phylotstratum ps
    with maximum e_c or e_ps is 1, and the relative
    expression levels of all other cells c or
    phylostrata ps range between 0 and 1.

    :param adata: The annotated data matrix of shape n_obs × n_vars. Rows correspond to cells and columns to genes.
    :param gene_id: Expects GeneID column from orthomap DataFrame.
    :param gene_age: Expects Phylostratum column from orthomap DataFrame.
    :param keep: Either define 'min' (ascending pre-sorting) or 'max' (non-ascending pre-sorting) to keep duplicates.
    :param layer: Layer to work on instead of X. If None, X is used.
    :param use: Specify if counts from adata.X (default) should be combined per age group to calculate
                the relative expression or if the corresponding 'pmatrix' (partial TEI values, see get_pmatrix)
                or 'wmatrix' (gene aged weighted expression) should be used.
                If layer is not None adata.X refers to adata.layers[layer].
    :param var_type: Specify how values should be combined per variable group. Possible values are 'mean', 'median',
                     'sum', 'min' and 'max'.
    :param group_by_obs: AnnData observation to be used as a group to combine count values.
    :param obs_fillna: Specify how NaN values should be named for observation.
    :param obs_type: Specify how values should be combined per observation group. Possible values are 'mean', 'median',
                     'sum', 'min' and 'max'.
    :param standard_scale: Wether or not to standardize the given axis (0: colums, 1: rows) between 0 and 1,
                           meaning for each variable or group, subtract the minimum and divide each by its maximum.
    :param normalize_total: Normalize counts per cell prior TEI calculation.
    :param log1p: Logarithmize the data matrix prior TEI calculation.
    :param target_sum: After normalization, each observation (cell) has a total count equal to target_sum.
    :return: Relative expression profile DataFrame.

    :type adata: AnnData
    :type gene_id: list
    :type gene_age: list
    :type keep: str
    :type layer: str
    :type use: str
    :type var_type: str
    :type group_by_obs: str
    :type obs_fillna: str
    :type obs_type: str
    :type standard_scale: int
    :type normalize_total: bool
    :type log1p: bool
    :type target_sum: float
    :rtype: pandas.DataFrame

    Example
    -------
    >>> import scanpy as sc
    >>> import matplotlib.pyplot as plt
    >>> import seaborn as sns
    >>> from orthomap import orthomap2tei, datasets
    >>> # download pre-calculated orthomap
    >>> #query_orthomap = orthomap2tei.read_orthomap(orthomapfile='Sun2021_Orthomap.tsv')
    >>> sun21_orthomap_file = datasets.sun21_orthomap(datapath='.')
    >>> # load query species orthomap
    >>> query_orthomap = orthomap2tei.read_orthomap(orthomapfile=sun21_orthomap_file)
    >>> # download and load scRNA data
    >>> #packer19_small = sc.read('packer19_small.h5ad')
    >>> packer19_small = datasets.packer19_small(datapath='.')
    >>> # get rematrix
    >>> packer19_small_rematrix = orthomap2tei.get_rematrix(
    >>>     adata=packer19_small,
    >>>     gene_id=query_orthomap['GeneID'],
    >>>     gene_age=query_orthomap['Phylostratum'])
    >>> # group by embryo.time.bin observation
    >>> packer19_small_rematrix_grouped = orthomap2tei.get_rematrix(
    >>>     adata=packer19_small,
    >>>     gene_id=query_orthomap['GeneID'],
    >>>     gene_age=query_orthomap['Phylostratum'],
    >>>     group_by_obs='embryo.time.bin')
    >>> # plot heatmap using partial TEI values
    >>> sns.heatmap(packer19_small_rematrix_grouped, cmap='viridis')
    >>> plt.show()
    >>> # group by embryo.time.bin observation and scale over rows
    >>> packer19_small_rematrix_grouped_rows = orthomap2tei.get_rematrix(
    >>>     adata=packer19_small,
    >>>     gene_id=query_orthomap['GeneID'],
    >>>     gene_age=query_orthomap['Phylostratum'],
    >>>     group_by_obs='embryo.time.bin',
    >>>     standard_scale=0)
    >>> # plot heatmap using partial TEI values
    >>> sns.heatmap(packer19_small_rematrix_grouped_rows, cmap='viridis')
    >>> plt.show()
    >>> # group by embryo.time.bin observation and scale over columns
    >>> packer19_small_rematrix_grouped_columns = orthomap2tei.get_rematrix(
    >>>     adata=packer19_small,
    >>>     gene_id=query_orthomap['GeneID'],
    >>>     gene_age=query_orthomap['Phylostratum'],
    >>>     group_by_obs='embryo.time.bin',
    >>>     standard_scale=1)
    >>> # plot heatmap using partial TEI values
    >>> sns.heatmap(packer19_small_rematrix_grouped_columns, cmap='viridis')
    >>> plt.show()
    """
    var_names_df,\
        id_age_df_keep_subset,\
        adata_counts,\
        var_names_subset,\
        sumx,\
        sumx_recd,\
        ps,\
        psd = _get_psd(adata=adata,
                       gene_id=gene_id,
                       gene_age=gene_age,
                       keep=keep,
                       layer=layer,
                       normalize_total=normalize_total,
                       log1p=log1p,
                       target_sum=target_sum)
    wmatrix = psd.dot(adata_counts.transpose()).transpose()
    pmatrix = sumx_recd.dot(wmatrix)
    tei = pmatrix.sum(1)
    phylostrata = list(set(id_age_df_keep_subset['Phylostrata']))
    rematrix = np.zeros((len(phylostrata), adata_counts.shape[0]))
    if use == 'pmatrix':
        for pk_idx, pk in enumerate(phylostrata):
            if var_type == 'mean':
                rematrix[pk_idx, ] = np.array(pmatrix[:, id_age_df_keep_subset['Phylostrata'].isin([pk])]
                                              .mean(1)).flatten()
            if var_type == 'median':
                rematrix[pk_idx, ] = np.apply_along_axis(
                    np.median, 1, pmatrix[:, id_age_df_keep_subset['Phylostrata'].isin([pk])].toarray()).flatten()
            if var_type == 'sum':
                rematrix[pk_idx, ] = np.array(pmatrix[:, id_age_df_keep_subset['Phylostrata'].isin([pk])]
                                              .sum(1)).flatten()
            if var_type == 'min':
                rematrix[pk_idx, ] = np.array(pmatrix[:, id_age_df_keep_subset['Phylostrata'].isin([pk])]
                                              .min(1).toarray()).flatten()
            if var_type == 'max':
                rematrix[pk_idx, ] = np.array(pmatrix[:, id_age_df_keep_subset['Phylostrata'].isin([pk])]
                                              .max(1).toarray()).flatten()
    elif use == 'wmatrix':
        for pk_idx, pk in enumerate(phylostrata):
            if var_type == 'mean':
                rematrix[pk_idx, ] = np.array(wmatrix[:, id_age_df_keep_subset['Phylostrata'].isin([pk])]
                                              .mean(1)).flatten()
            if var_type == 'median':
                rematrix[pk_idx, ] = np.apply_along_axis(
                    np.median, 1, wmatrix[:, id_age_df_keep_subset['Phylostrata'].isin([pk])].toarray()).flatten()
            if var_type == 'sum':
                rematrix[pk_idx, ] = np.array(wmatrix[:, id_age_df_keep_subset['Phylostrata'].isin([pk])]
                                              .sum(1)).flatten()
            if var_type == 'min':
                rematrix[pk_idx, ] = np.array(wmatrix[:, id_age_df_keep_subset['Phylostrata'].isin([pk])]
                                              .min(1).toarray()).flatten()
            if var_type == 'max':
                rematrix[pk_idx, ] = np.array(wmatrix[:, id_age_df_keep_subset['Phylostrata'].isin([pk])]
                                              .max(1).toarray()).flatten()
    else:
        for pk_idx, pk in enumerate(phylostrata):
            if var_type == 'mean':
                rematrix[pk_idx, ] = np.array(adata_counts[:, id_age_df_keep_subset['Phylostrata'].isin([pk])]
                                              .mean(1)).flatten()
            if var_type == 'median':
                rematrix[pk_idx, ] = np.apply_along_axis(
                    np.median, 1, adata_counts[:, id_age_df_keep_subset['Phylostrata'].isin([pk])].toarray()).flatten()
            if var_type == 'sum':
                rematrix[pk_idx, ] = np.array(adata_counts[:, id_age_df_keep_subset['Phylostrata'].isin([pk])]
                                              .sum(1)).flatten()
            if var_type == 'min':
                rematrix[pk_idx, ] = np.array(adata_counts[:, id_age_df_keep_subset['Phylostrata'].isin([pk])]
                                              .min(1).toarray()).flatten()
            if var_type == 'max':
                rematrix[pk_idx, ] = np.array(adata_counts[:, id_age_df_keep_subset['Phylostrata'].isin([pk])]
                                              .max(1).toarray()).flatten()
    rematrix_df = pd.DataFrame(rematrix)
    rematrix_df['ps'] = phylostrata
    rematrix_df.set_index('ps',
                          inplace=True)
    rematrix_df.columns = adata.obs_names
    if group_by_obs is not None:
        if adata.obs[group_by_obs].dtype.name == 'category':
            obs_group_nan = pd.DataFrame(adata.obs[group_by_obs])
        else:
            obs_group_nan = pd.DataFrame(adata.obs[group_by_obs].fillna(obs_fillna))
        if obs_type == 'mean':
            rematrix_df = \
                rematrix_df.transpose().groupby(obs_group_nan[group_by_obs]).mean().transpose()
        if obs_type == 'median':
            rematrix_df = \
                rematrix_df.transpose().groupby(obs_group_nan[group_by_obs]).median().transpose()
        if obs_type == 'sum':
            rematrix_df = \
                rematrix_df.transpose().groupby(obs_group_nan[group_by_obs]).sum().transpose()
        if obs_type == 'min':
            rematrix_df = \
                rematrix_df.transpose().groupby(obs_group_nan[group_by_obs]).min().transpose()
        if obs_type == 'max':
            rematrix_df = \
                rematrix_df.transpose().groupby(obs_group_nan[group_by_obs]).max().transpose()
    if standard_scale is not None:
        if standard_scale == 0:
            rematrix_df = rematrix_df.apply(_min_max_to_01,
                                            axis=1,
                                            raw=True)
        if standard_scale == 1:
            rematrix_df = rematrix_df.apply(_min_max_to_01,
                                            axis=0,
                                            raw=True)
    return rematrix_df


def _get_min_max_array(ndarray,
                       min_expr=None,
                       max_expr=None):
    """
    A helper function to reduce data between min_expr and max_expr.

    :param ndarray: Data to be reduced.
    :param min_expr: Specify minimal expression to be included.
    :param max_expr: Specify maximal expression to be included.
    :return: Reduced data.

    :type ndarray:
    :type min_expr: float
    :type max_expr: float
    :rtype: numpy.ndarray

    Example
    -------
    >>> import numpy as np
    >>> from orthomap import orthomap2tei
    >>> data = np.array([0, 0, 1, 1, 4, 3, 7, 0, 10, 12])
    >>> orthomap2tei._get_min_max_array(ndarray=data, min_expr=1, max_expr=10)
    """
    if min_expr is not None and max_expr is not None:
        return ndarray[np.where(np.logical_and(ndarray >= min_expr,
                                               ndarray <= max_expr))]
    if min_expr is not None and max_expr is None:
        return ndarray[np.where(ndarray >= min_expr)]
    if min_expr is None and max_expr is not None:
        return ndarray[np.where(ndarray <= max_expr)]
    else:
        return ndarray


def mergeby_from_dataframe(df,
                           col_group=None,
                           col_fillna='__NaN',
                           row_group=None,
                           row_fillna='__NaN',
                           level='col',
                           min_expr=None,
                           max_expr=None):
    """
    This function groups all values of an existing pandas.DataFrame as an array based on column or row groups.
    The resulting pandas.DataFrame can be used to e.g. apply statistics or visualize the groups more easily.

    :param df: DataFrame to be grouped
    :param col_group: numpy.ndarray containing group assignment for columns
                      (needs to be same length and same order as Columns of df).
    :param col_fillna:
    :param row_group: numpy.ndarray containing group assignment for row
                      (needs to be same length and same order as Index of df).
    :param row_fillna:
    :param level: Specify if col or row should be used as primary group (only effects output orientation).
    :param min_expr: Specify minimal expression to be included.
    :param max_expr: Specify maximal expression to be included.
    :return: List of three DataFrame. First DataFrame contains the grouped data (each cell contains a numpy.ndarray).
             Second DataFrame contains original col and third row assignment and groupings.

    :type df: pandas.DataFrame
    :type col_group: numpy.ndarray
    :type col_fillna: str
    :type row_group: numpy.ndarray
    :type row_fillna: str
    :type level: str
    :type min_expr: float
    :type max_expr: float
    :rtype: list
    """
    merged_df = None
    col_assignment = pd.DataFrame(list(df.columns), columns=['orig.col'])
    row_assignment = pd.DataFrame(list(df.index), columns=['orig.row'])
    if col_group is not None:
        col_group_df = pd.DataFrame(col_group).fillna(col_fillna)
        col_group_df.columns = ['col_group']
        if col_group_df.shape[0] != df.shape[1]:
            print('col_group length and number of columns of df differ')
            return
        col_grouped = col_group_df.groupby('col_group')
        col_groups = col_grouped.groups.keys()
        col_assignment['col_group'] = col_group_df['col_group']
    if row_group is not None:
        row_group_df = pd.DataFrame(row_group).fillna(row_fillna)
        row_group_df.columns = ['row_group']
        if row_group_df.shape[0] != df.shape[0]:
            print('row_group length and number of rows of df differ')
            return
        row_grouped = row_group_df.groupby('row_group')
        row_groups = row_grouped.groups.keys()
        row_assignment['row_group'] = row_group_df['row_group']
    if col_group is not None and row_group is not None:
        merged_df = pd.DataFrame(row_groups).set_index(0)
        for col_idx, (col_grp, col_grp_idx) in enumerate(col_grouped.indices.items()):
            merged_df[col_grp] = [[] for x in row_groups]
            for row_idx, (row_grp, row_grp_idx) in enumerate(row_grouped.indices.items()):
                merged_df.iat[row_idx, col_idx] = _get_min_max_array(df.iloc[row_grp_idx,
                                                                             col_grp_idx].values.flatten(),
                                                                     min_expr=min_expr,
                                                                     max_expr=max_expr)
        if level == 'col':
            return [merged_df,
                    col_assignment,
                    row_assignment]
        if level == 'row':
            return [merged_df.transpose(),
                    col_assignment,
                    row_assignment]
    if col_group is None and row_group is not None:
        merged_df = pd.DataFrame(df.columns.values).set_index(0)
        for row_idx, (row_grp, row_grp_idx) in enumerate(row_grouped.indices.items()):
            merged_df[row_grp] = df.transpose().iloc[:, row_grp_idx].apply(np.array, axis=1).apply(_get_min_max_array,
                                                                                                   min_expr=min_expr,
                                                                                                   max_expr=max_expr)
        if level == 'col':
            return [merged_df.transpose(),
                    col_assignment,
                    row_assignment]
        if level == 'row':
            return [merged_df,
                    col_assignment,
                    row_assignment]
    if col_group is not None and row_group is None:
        merged_df = pd.DataFrame(df.index).set_index(df.index.name)
        for col_idx, (col_grp, col_grp_idx) in enumerate(col_grouped.indices.items()):
            merged_df[col_grp] = df.iloc[:, col_grp_idx].apply(np.array, axis=1).apply(_get_min_max_array,
                                                                                       min_expr=min_expr,
                                                                                       max_expr=max_expr)
        if level == 'col':
            return [merged_df,
                    col_assignment,
                    row_assignment]
        if level == 'row':
            return [merged_df.transpose(),
                    col_assignment,
                    row_assignment]
    if level == 'col':
        return [df,
                col_assignment,
                row_assignment]
    if level == 'row':
        return [df.transpose(),
                col_assignment,
                row_assignment]


def mergeby_from_counts(adata,
                        layer=None,
                        group_by_var=None,
                        var_fillna='__NaN',
                        group_by_obs=None,
                        obs_fillna='__NaN',
                        level='obs',
                        min_expr=None,
                        max_expr=None,
                        normalize_total=False,
                        log1p=False,
                        target_sum=1e6):
    """
    This function groups all counts of an existing AnnData object as an array based on variable or observation groups.
    The resulting pandas.DataFrame can be used to e.g. apply statistics or visualize the groups more easily.

    :param adata: The annotated data matrix of shape n_obs × n_vars. Rows correspond to cells and columns to genes.
    :param layer: Layer to work on instead of X. If None, X is used.
    :param group_by_var: AnnData variable to be used as a group to combine count values.
    :param var_fillna: Specify how NaN values should be named for variable.
    :param group_by_obs: AnnData observation to be used as a group to combine count values.
    :param obs_fillna: Specify how NaN values should be named for observation.
    :param level: Specify if observation or variable should be used as primary group (only effects output orientation).
    :param min_expr: Specify minimal expression to be included.
    :param max_expr: Specify maximal expression to be included.
    :param normalize_total: Normalize counts per cell.
    :param log1p: Logarithmize the data matrix.
    :param target_sum: After normalization, each observation (cell) has a total count equal to target_sum.
    :return: List of three DataFrame. First DataFrame contains the grouped data (each cell contains a numpy.ndarray).
             Second DataFrame contains original variable and observation assignment and groupings.

    :type adata: AnnData
    :type layer: str
    :type group_by_var: str
    :type var_fillna: str
    :type group_by_obs: str
    :type obs_fillna: str
    :type level: str
    :type min_expr: float
    :type max_expr: float
    :type normalize_total: bool
    :type log1p: bool
    :type target_sum: float
    :rtype: list

    Example
    -------
    >>> import scanpy as sc
    >>> import matplotlib.pyplot as plt
    >>> import seaborn as sns
    >>> from orthomap import orthomap2tei, datasets
    >>> # download and load scRNA data
    >>> #packer19_small = sc.read('packer19_small.h5ad')
    >>> packer19_small = datasets.packer19_small(datapath='.')
    >>> # add gene age values to existing adata object
    >>> orthomap2tei.add_gene_age2adata_var(
    >>>     adata=packer19_small,
    >>>     gene_id=query_orthomap['GeneID'],
    >>>     gene_age=query_orthomap['Phylostratum'])
    >>> # get group counts
    >>> packer19_small_group_counts = orthomap2tei.mergeby_from_counts(
    >>>     adata=packer19_small,
    >>>     group_by_var='Phylostrata',
    >>>     group_by_obs='embryo.time.bin')
    """
    adata_counts = _get_counts(adata=adata,
                               layer=layer,
                               normalize_total=normalize_total,
                               log1p=log1p,
                               target_sum=target_sum)
    merged_df = None
    var_assignment = pd.DataFrame(list(adata.var.index), columns=['orig.var'])
    obs_assignment = pd.DataFrame(list(adata.obs.index), columns=['orig.obs'])
    if group_by_obs is not None:
        if adata.obs[group_by_obs].dtype.name == 'category':
            obs_group_df = pd.DataFrame(adata.obs[group_by_obs])
        else:
            obs_group_df = pd.DataFrame(adata.obs[group_by_obs]).fillna(obs_fillna)
        if obs_group_df.shape[0] != adata_counts.shape[0]:
            print('obs_group length and number of cells of adata_counts differ')
            return
        obs_grouped = obs_group_df.groupby(group_by_obs)
        obs_groups = obs_grouped.groups.keys()
        obs_assignment[group_by_obs] = np.array(obs_group_df[group_by_obs])
    if group_by_var is not None:
        if adata.obs[group_by_var].dtype.name == 'category':
            var_group_df = pd.DataFrame(adata.obs[group_by_var])
        else:
            var_group_df = pd.DataFrame(adata.var[group_by_var]).fillna(var_fillna)
        if var_group_df.shape[0] != adata_counts.shape[1]:
            print('var_group length and number of genes of adata_counts differ')
            return
        var_grouped = var_group_df.groupby(group_by_var)
        var_groups = var_grouped.groups.keys()
        var_assignment[group_by_var] = np.array(var_group_df[group_by_var])
    if group_by_obs is not None and group_by_var is not None:
        merged_df = pd.DataFrame(var_groups).set_index(0)
        for obs_idx, (obs_grp, obs_grp_idx) in enumerate(obs_grouped.indices.items()):
            merged_df[obs_grp] = [[] for x in var_groups]
            for var_idx, (var_grp, var_grp_idx) in enumerate(var_grouped.indices.items()):
                group_counts = _get_min_max_array(ndarray=np.array(adata_counts[:, var_grp_idx]
                                                                   [obs_grp_idx, :].toarray()).flatten(),
                                                  min_expr=min_expr,
                                                  max_expr=max_expr)
                merged_df.iat[var_idx, obs_idx] = group_counts
        if level == 'obs':
            return [merged_df.transpose(),
                    obs_assignment,
                    var_assignment]
        if level == 'row':
            return [merged_df,
                    obs_assignment,
                    var_assignment]
    if group_by_obs is not None and group_by_var is None:
        merged_df = pd.DataFrame(list(obs_groups)).set_index(0)
        for obs_idx, (obs_grp, obs_grp_idx) in enumerate(obs_grouped.indices.items()):
            group_counts = _get_min_max_array(ndarray=np.array(adata_counts[obs_grp_idx, :].toarray()).flatten(),
                                              min_expr=min_expr,
                                              max_expr=max_expr)
            merged_df[obs_grp] = group_counts
        if level == 'obs':
            return [merged_df.transpose(),
                    obs_assignment,
                    var_assignment]
        if level == 'var':
            return [merged_df,
                    obs_assignment,
                    var_assignment]
    if group_by_obs is None and group_by_var is not None:
        merged_df = pd.DataFrame(list(var_groups)).set_index(0)
        for var_idx, (var_grp, var_grp_idx) in enumerate(var_grouped.indices.items()):
            group_counts = _get_min_max_array(ndarray=np.array(adata_counts[:, var_grp_idx].toarray()).flatten(),
                                              min_expr=min_expr,
                                              max_expr=max_expr)
            merged_df[var_grp] = group_counts
        if level == 'obs':
            return [merged_df.transpose(),
                    obs_assignment,
                    var_assignment]
        if level == 'var':
            return [merged_df,
                    obs_assignment,
                    var_assignment]
    if level == 'obs':
        return [pd.DataFrame.sparse.from_spmatrix(adata_counts).transpose(),
                obs_assignment,
                var_assignment]
    if level == 'var':
        return [pd.DataFrame.sparse.from_spmatrix(adata_counts),
                obs_assignment,
                var_assignment]


def _get_min_max_expr_number(ndarray,
                             min_expr=1,
                             max_expr=None):
    """
    A helper function to

    :param ndarray:
    :param min_expr:
    :param max_expr:
    :return:

    Example
    --------
    >>>
    """
    if type(min_expr) == str:
        min_expr = np.quantile(a=ndarray, q=float(min_expr.split('q')[1])/100)
    if type(max_expr) == str:
        max_expr = np.quantile(a=ndarray, q=float(max_expr.split('q')[1])/100)
    if max_expr is not None:
        return np.bitwise_and(np.greater_equal(ndarray, min_expr), np.less(ndarray, max_expr)).sum()
    else:
        return np.greater_equal(ndarray, min_expr).sum()


def _get_quantile_expr_number(ndarray,
                              quantile=[0, 5, 25, 50, 75, 95],
                              min_expr=1,
                              max_expr=None):
    """
    A helper function to

    :param ndarray:
    :param quantile:
    :param min_expr:
    :param max_expr:
    :return:

    Example
    --------
    >>>
    """
    q_expr = []
    if type(min_expr) == str:
        min_expr = np.quantile(a=ndarray, q=float(min_expr.split('q')[1])/100)
    if type(max_expr) == str:
        max_expr = np.quantile(a=ndarray, q=float(max_expr.split('q')[1])/100)
    for q in quantile:
        if max_expr is not None:
            q_expr.append(np.percentile(a=ndarray[np.bitwise_and(np.greater_equal(ndarray, min_expr),
                                                                 np.less(ndarray, max_expr))], q=q))
        else:
            q_expr.append(np.percentile(a=ndarray[np.greater_equal(ndarray, min_expr)], q=q))
    q_expr_number = []
    for qe in q_expr:
        q_expr_number.append(_get_min_max_expr_number(ndarray, min_expr=qe, max_expr=max_expr))
    return q_expr_number


def get_e50(adata,
            gene_id,
            gene_age,
            keep='min',
            layer=None,
            group_by_var=None,
            var_type='mean',
            group_by_obs=None,
            obs_type='mean',
            standard_scale=None,
            normalize_total=False,
            log1p=False,
            target_sum=1e6,
            min_expr=1,
            max_expr=None):
    """

    :param adata:
    :param gene_id:
    :param gene_age:
    :param keep:
    :param layer:
    :param use:
    :param col_type:
    :param standard_scale:
    :param group_by:
    :param group_type:
    :param normalize_total:
    :param log1p:
    :param target_sum:
    :param min_expr:
    :param max_expr:
    :return:

    Example
    --------
    >>>
    """
    #id_age_df_keep_subset, adata_counts, var_names_subset, sumx, sumx_recd, ps, psd =\
    #    _get_psd(adata, gene_id, gene_age, keep, layer, normalize_total, log1p, target_sum)
    #phylostrata = list(set(id_age_df_keep_subset['Phylostrata']))
    #min_expr_global = np.apply_along_axis(
    #    _get_min_max_expr_number, 1, adata_counts.toarray(), min_expr=min_expr)
    #e50_global = np.apply_along_axis(
    #    _get_quantile_expr_number, 1, adata_counts.toarray(), quantile=50, min_expr=min_expr)
    #for pk_idx, pk in enumerate(phylostrata):
    #    e50matrix[pk_idx,] = np.apply_along_axis(
    #        np.median, 1, adata_counts[:, id_age_df_keep_subset['Phylostrata'].isin([pk])].toarray()).flatten()
    #if group_by is not None:
    #if standard_scale is not None:
    return


def get_bins(tobin_df,
             bincol,
             q=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9],
             method='median_unbiased'):
    """
    This function sorts values of a DataFrame column and return the binary categories in which they fall.

    :param tobin_df: DataFrame which contains the column that should be binned.
    :param bincol: Name of the columns that should be binned.
    :param q: Array of inner quantiles to be used for binning. Border quantiles will be set automatically.
    :param method: This parameter specifies the method to use for estimating the quantile.
    :return: DataFrame with additional column which contain the binned values as categories.

    :type tobin_df: pandas.DataFrame
    :type bincol: str
    :type q: list
    :type method: str
    :rytpe: pandas.DataFrame
    """
    if len(q) == 0:
        print('at least provide one quantile value')
    tobin_df[bincol + '_binned'] = np.nan
    tobin_df[bincol + '_bins'] = np.nan
    qs = np.nanquantile(a=tobin_df[bincol], q=q, method=method)
    for qs_idx, qs_val in enumerate(qs):
        if qs_idx == 0:
            if qs[qs_idx] == qs[qs_idx + 1]:
                print('bin cutoff neighbors are equal, please reduce bin size')
            tobin_df.loc[tobin_df[bincol] < qs[qs_idx], [bincol + '_binned']] = qs_idx + 1
            tobin_df.loc[tobin_df[bincol] < qs[qs_idx], [bincol + '_bins']] = 'x < ' + str(qs[qs_idx])
            tobin_df.loc[np.bitwise_and(tobin_df[bincol] >= qs[qs_idx], tobin_df[bincol] < qs[qs_idx + 1]), [
                bincol + '_binned']] = qs_idx + 2
            tobin_df.loc[np.bitwise_and(tobin_df[bincol] >= qs[qs_idx], tobin_df[bincol] < qs[qs_idx + 1]), [
                bincol + '_bins']] = str(qs[qs_idx]) + ' >= x < ' + str(qs[qs_idx + 1])
        elif qs_idx+1 == len(qs):
            tobin_df.loc[tobin_df[bincol] >= qs[qs_idx], [bincol + '_binned']] = qs_idx + 2
            tobin_df.loc[tobin_df[bincol] >= qs[qs_idx], [bincol + '_bins']] = str(qs[qs_idx]) + ' < x'
        else:
            if qs[qs_idx] == qs[qs_idx + 1]:
                print('bin cutoff neighbors are equal, please reduce bin size')
            tobin_df.loc[np.bitwise_and(tobin_df[bincol] >= qs[qs_idx], tobin_df[bincol] < qs[qs_idx + 1]), [
                bincol + '_binned']] = qs_idx + 2
            tobin_df.loc[np.bitwise_and(tobin_df[bincol] >= qs[qs_idx], tobin_df[bincol] < qs[qs_idx + 1]), [
                bincol + '_bins']] = str(qs[qs_idx]) + ' >= x < ' + str(qs[qs_idx + 1])
    return tobin_df
