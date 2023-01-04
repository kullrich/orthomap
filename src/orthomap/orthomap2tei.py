"""
Author: Kristian K Ullrich
date: November 2022
email: ullrich@evolbio.mpg.de
License: GPL-3
"""


import os

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy
from alive_progress import alive_bar


def read_orthomap(orthomapfile):
    """

    :param orthomapfile: str
        File name of pre-calculated orthomap file.
    :return:

    Example
    --------
    >>> from orthomap import orthomap2tei
    >>> # load query species orthomap
    >>> query_orthomap = orthomap2tei.read_orthomap('Sun2021_Orthomap.tsv')
    >>> query_orthomap
    """
    if os.path.exists(orthomapfile):
        orthomap = pd.read_csv(orthomapfile, delimiter='\t')
    return orthomap


def geneset_overlap(geneset1, geneset2):
    """

    :param geneset1: list
        List of gene or transcript names set 1.
    :param geneset2: list
        List of gene or transcript names set 2.
    :return:

    Example
    --------
    >>> from orthomap import orthomap2tei
    >>> g1 = ['g1.1', 'g1.2', 'g2.1', 'g3.1', 'g3.2']
    >>> g2 = ['g1.1', 'g2.1', 'g3.1']
    >>> orthomap2tei.geneset_overlap(g1, g2)
    """
    g1 = set(geneset1)
    g2 = set(geneset2)
    g1_g2 = g1.intersection(g2)
    d = {'g1_g2_overlap': [len(g1_g2)],
         'g1_ratio': [len(g1_g2)/len(g1)],
         'g2_ratio': [len(g1_g2)/len(g2)]}
    df = pd.DataFrame(data=d)
    return df


def replace_by(x_orig, xmatch, xreplace):
    """

    :param x_orig:
    :param xmatch:
    :param xreplace:
    :return:

    Example
    --------
    >>>
    """
    replace_dict = {}
    for i, j in enumerate(xmatch):
        replace_dict[j] = xreplace[i]
    x_new = [replace_dict[x] for x in x_orig]
    return x_new


def keep_min_max(df, keep='min', dup_col=['GeneID'], sort_col=['Phylostrata']):
    """

    :param df: DataFrame
        Expects orthomap DataFrame, but can be any if dup_col and sort_col are present.
    :param keep: str (default: min)
        Either define 'min' (ascending pre-sorting) or 'max' (non-ascending pre-sorting) to keep duplicates.
    :param dup_col: [str] (default: GeneID)
        Column name(s) to be searched for duplicates.
    :param sort_col: [str] (default: Phylostrata)
        Column names(s) to sort DataFrame.
    :return:

    Example
    --------
    >>> from orthomap import orthomap2tei
    >>> # create artificial DataFrame
    >>> import pandas as pd
    >>> my_orthomap = pd.DataFrame.from_dict({'GeneID':['g1', 'g1', 'g2', 'g3', 'g3'],\
    >>> 'Phylostrata':[3, 1, 2, 5, 7]})
    >>> # keep min value
    >>> orthomap2tei.keep_min_max(my_orthomap, keep='min')
    >>> # keep max value
    >>> orthomap2tei.keep_min_max(my_orthomap, keep='max')
    """
    if keep == 'min':
        df_sorted = df.sort_values(by=sort_col, ascending=True)
    if keep == 'max':
        df_sorted = df.sort_values(by=sort_col, ascending=False)
    keep_idx = ~df_sorted.duplicated(dup_col)
    df_out = df_sorted[keep_idx]
    return df_out


def split_gene_id_by_gene_age(gene_id, gene_age, keep='min', adata=None):
    """

    :param gene_id: list
        Expects GeneID column from orthomap DataFrame.
    :param gene_age: list
        Expects Phylostratum column from orthomap DataFrame.
    :param keep: str (default: min)
        In case of duplicated GeneIDs with different Phylostrata assignments, either keep 'min' or 'max' value.
    :param adata: AnnData (default: None)
        The annotated data matrix of shape n_obs × n_vars. Rows correspond to cells and columns to genes.
    :return:
    """
    id_age_df = pd.DataFrame(data={'GeneID': gene_id, 'Phylostrata': gene_age})
    # check and drop duplicated GeneID
    id_age_df_keep = keep_min_max(id_age_df, keep=keep, dup_col=['GeneID'], sort_col=['Phylostrata'])
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
        gene_id_gene_age_dict[str(pk)] = list(id_age_df_keep_subset[id_age_df_keep_subset['Phylostrata'] == pk]['GeneID'])
    return gene_id_gene_age_dict


def get_psd(adata, gene_id, gene_age, keep='min', layer=None, normalize_total=False, log1p=False, target_sum=1e6):
    """

    :param adata: AnnData
        The annotated data matrix of shape n_obs × n_vars. Rows correspond to cells and columns to genes.
    :param gene_id: list
        Expects GeneID column from orthomap DataFrame.
    :param gene_age:
        Expects Phylostratum column from orthomap DataFrame.
    :param keep: str (default: min)
        In case of duplicated GeneIDs with different Phylostrata assignments, either keep 'min' or 'max' value.
    :param layer: Optional[str] (default: None)
        Layer to work on instead of X. If None, X is used.
    :param normalize_total: bool (default: False)
        Normalize counts per cell prior TEI calculation.
    :param log1p: bool (default: False)
        Logarithmize the data matrix prior TEI calculation.
    :param target_sum: Optional[float] (default: 1e6)
        After normalization, each observation (cell) has a total count equal to target_sum.
    :return:

    Example
    --------
    >>> from orthomap import orthomap2tei
    >>> # load query species orthomap
    >>> query_orthomap = orthomap2tei.read_orthomap('Sun2021_Orthomap.tsv')
    >>> # load scRNA data
    >>> import scanpy as sc
    >>> celegans_data = sc.read('celegans.h5ad')
    >>> # get psd from existing adata object
    >>> celegans_id_age_df_keep_subset, celegans_adata_counts, celegans_var_names_subset, celegans_sumx,\
    >>> celegans_sumx_recd, celegans_ps, celegans_psd =\
    >>> orthomap2tei.get_psd(adata=celegans_data,\
    >>> gene_id=query_orthomap['GeneID'],\
    >>> gene_age=query_orthomap['Phylostratum'])
    """
    id_age_df = pd.DataFrame(data={'GeneID': gene_id, 'Phylostrata': gene_age})
    # check and drop duplicated GeneID
    id_age_df_keep = keep_min_max(id_age_df, keep=keep, dup_col=['GeneID'], sort_col=['Phylostrata'])
    # get overlap
    gene_intersection = pd.Index(id_age_df_keep['GeneID']).intersection(adata.var_names)
    id_age_df_keep_subset = id_age_df_keep.loc[id_age_df_keep['GeneID'].isin(gene_intersection)]
    id_age_df_keep_subset = id_age_df_keep_subset.sort_values('GeneID')
    adata_counts = adata.X
    if layer is not None:
        adata_counts = adata.layers[layer]
    if normalize_total and log1p:
        adata_norm = sc.pp.normalize_total(adata, target_sum=target_sum, layer=layer, copy=True)
        sc.pp.log1p(adata_norm, layer=layer)
        adata_counts = adata_norm.X
        if layer is not None:
            adata_counts = adata_norm.layers[layer]
    if normalize_total and not log1p:
        adata_norm = sc.pp.normalize_total(adata, target_sum=target_sum, layer=layer, copy=True)
        adata_counts = adata_norm.X
        if layer is not None:
            adata_counts = adata_norm.layers[layer]
    if not normalize_total and log1p:
        adata_log1p = sc.pp.log1p(adata, layer=layer, copy=True)
        adata_counts = adata_log1p.X
        if layer is not None:
            adata_counts = adata_log1p.layers[layer]
    adata_counts = adata_counts[:, adata.var_names.isin(id_age_df_keep_subset['GeneID'])]
    var_names_subset = adata.var_names[adata.var_names.isin(id_age_df_keep_subset['GeneID'])]
    var_names_subset_idx = var_names_subset.sort_values(return_indexer=True)[1]
    adata_counts = adata_counts[:, var_names_subset_idx]
    sumx = adata_counts.sum(1)
    sumx_rec = np.reciprocal(sumx)
    sumx_recd = scipy.sparse.diags(np.array(sumx_rec).flatten())
    ps = np.array(id_age_df_keep_subset['Phylostrata'])
    psd = scipy.sparse.diags(ps)
    return [id_age_df_keep_subset, adata_counts, var_names_subset, sumx, sumx_recd, ps, psd]


def get_tei(adata, gene_id, gene_age, keep='min', layer=None, add=True, obs_name='tei', boot=False, bt=10,
            normalize_total=False, log1p=False, target_sum=1e6):
    """
    This function computes the phylogenetically based transcriptome evolutionary
    index (TEI) similar to Domazet-Loso & Tautz, 2010.

     The TEI measure represents the weighted arithmetic mean
     (expression levels as weights for the phylostratum value) over all
     evolutionary age categories denoted as _phylostra_.

     \deqn{TEI_s = sum (e_is * ps_i) / sum e_is}

     where \eqn{TEI_s} denotes the TEI value in developmental stage \eqn{s, e_is}
     denotes the gene expression level of gene \eqn{i} in stage \eqn{s}, and \eqn{ps_i}
     denotes the corresponding phylostratum of gene \eqn{i, i = 1,...,N} and
     \eqn{N = total number of genes}.

    If the parameter boot is set to true,
    the strata values are sampled and the global TEI
    is calculated bt times.

    :param adata: AnnData
        The annotated data matrix of shape n_obs × n_vars. Rows correspond to cells and columns to genes.
    :param gene_id: list
        Expects GeneID column from orthomap DataFrame.
    :param gene_age: list
        Expects GeneID column from orthomap DataFrame.
    :param keep: str (default: min)
        In case of duplicated GeneIDs with different Phylostrata assignments, either keep 'min' or 'max' value.
    :param layer: Optional[str] (default: None)
        Layer to work on instead of X. If None, X is used.
    :param add: bool (default: True)
        Add TEI values as observation to existing adata object using obs_name.
    :param obs_name: str
        Observation name to be used for TEI values in existing adata object.
    :param boot: bool (default: False)
       Specify if bootstrap TEI values should be calculated and returned as DataFrame.
    :param bt: int (default: 10)
        Number of bootstrap to calculate.
    :param normalize_total: bool (default: False)
        Normalize counts per cell prior TEI calculation.
    :param log1p: bool (default: False)
        Logarithmize the data matrix prior TEI calculation.
    :param target_sum: Optional[float] (default: 1e6)
        After normalization, each observation (cell) has a total count equal to target_sum.
    :return:

    Example
    --------
    >>> from orthomap import orthomap2tei
    >>> # load query species orthomap
    >>> query_orthomap = orthomap2tei.read_orthomap('Sun2021_Orthomap.tsv')
    >>> # load scRNA data
    >>> import scanpy as sc
    >>> celegans_data = sc.read('celegans.h5ad')
    >>> # add TEI values to existing adata object
    >>> orthomap2tei.get_tei(adata=celegans_data,\
    >>> gene_id=query_orthomap['GeneID'],\
    >>> gene_age=query_orthomap['Phylostratum'],\
    >>> add=True)
    >>> # get 10 bootstap TEI values
    >>> orthomap2tei.get_tei(adata=celegans_data,\
    >>> gene_id=query_orthomap['GeneID'],\
    >>> gene_age=query_orthomap['Phylostratum'],\
    >>> boot=True, bt=10)
    """
    id_age_df_keep_subset, adata_counts, var_names_subset, sumx, sumx_recd, ps, psd =\
        get_psd(adata, gene_id, gene_age, keep, layer, normalize_total, log1p, target_sum)
    teimatrix = psd.dot(adata_counts.transpose()).transpose()
    pmatrix = sumx_recd.dot(teimatrix)
    tei = pmatrix.sum(1)
    tei_df = pd.DataFrame(tei, columns=['tei'])
    tei_df.index = adata.obs_names
    if add:
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


def get_pmatrix(adata, gene_id, gene_age, keep='min', layer=None, layer_name='pmatrix', add_obs=True, add_var=True,
                normalize_total=False, log1p=False, target_sum=1e6):
    """
    This function computes the partial transcriptome evolutionary index (TEI) values for each single gene.

    In detail, each gene gets a TEI contribution profile as follows:
    \deqn{TEI_is = f_is * ps_i}
    where TEI_is is the partial TEI value of gene i,
    \eqn{f_is = e_is / \sum e_is} and \eqn{ps_i} is the phylostratum of gene i.

    The partial TEI matrix can be used to perform different cluster
    analyses and also gives an overall impression of the contribution of each
    gene to the global TEI pattern.

    :param adata: AnnData The annotated data matrix of shape n_obs × n_vars. Rows correspond to cells and columns to genes.
    :param gene_id: list
        Expects GeneID column from orthomap DataFrame.
    :param gene_age: list
        Expects Phylostratum column from orthomap DataFrame.
    :param keep: str (default: min)
        Either define 'min' (ascending pre-sorting) or 'max' (non-ascending pre-sorting) to keep duplicates.
    :param layer: Optional[str] (default: None)
        Layer to work on instead of X. If None, X is used.
    :param layer_name: str (default: pmatrix)
        Layer to add to existing adata object. If copy, return a new adata object.
    :param add_obs:
    :param add_var:
    :param normalize_total: bool (default: False)
        Normalize counts per cell prior TEI calculation.
    :param log1p: bool (default: False)
        Logarithmize the data matrix prior TEI calculation.
    :param target_sum: Optional[float] (default: 1e6)
        After normalization, each observation (cell) has a total count equal to target_sum.
    :return:

    Example
    --------
    >>>
    """
    id_age_df_keep_subset, adata_counts, var_names_subset, sumx, sumx_recd, ps, psd =\
        get_psd(adata, gene_id, gene_age, keep, layer, normalize_total, log1p, target_sum)
    teimatrix = psd.dot(adata_counts.transpose()).transpose()
    pmatrix = sumx_recd.dot(teimatrix)
    tei = pmatrix.sum(1)
    adata_pmatrix = ad.AnnData(adata_counts)
    adata_pmatrix.layers[layer_name] = pmatrix
    adata_pmatrix.obs_names = adata.obs_names
    adata_pmatrix.var_names = var_names_subset
    if add_var:
        for kv in adata.var.keys():
            adata_pmatrix.var[kv] = adata.var[kv][adata.var_names.isin(id_age_df_keep_subset['GeneID'])]
    if add_obs:
        for ko in adata.obs.keys():
            adata_pmatrix.obs[ko] = adata.obs[ko]
    return adata_pmatrix


def get_pstrata(adata, gene_id, gene_age, keep='min', layer=None, cumsum=False, standard_scale=None, group_by=None,
                normalize_total=False, log1p=False, target_sum=1e6):
    """
    This function computes the partial transcriptome evolutionary index (TEI) values combined for each strata.

    In detail, each gene gets a TEI contribution profile as follows:
    \deqn{TEI_is = f_is * ps_i}
    where TEI_is is the partial TEI value of gene i,
    \eqn{f_is = e_is / \sum e_is} and \eqn{ps_i} is the phylostratum of gene i.

    \eqn{TEI_is} values are combined per \eqn{ps}.

    The partial TEI values combined per strata give an overall impression
    of the contribution of each strata to the global TEI pattern.

    :param adata: AnnData The annotated data matrix of shape n_obs × n_vars. Rows correspond to cells and columns to genes.
    :param gene_id: list
        Expects GeneID column from orthomap DataFrame.
    :param gene_age: list
        Expects Phylostratum column from orthomap DataFrame.
    :param keep: str (default: min)
        Either define 'min' (ascending pre-sorting) or 'max' (non-ascending pre-sorting) to keep duplicates.
    :param layer: Optional[str] (default: None)
        Layer to work on instead of X. If None, X is used.
    :param cumsum:
    :param standard_scale: int (default: None)
        Whether or not to standardize the given axis (0: colums, 1: rows) between 0 and 1,
        meaning for each variable or group, subtract the minimum and divide each by its maximum.
    :param group_by:
    :param normalize_total: bool (default: False)
        Normalize counts per cell prior TEI calculation.
    :param log1p: bool (default: False)
        Logarithmize the data matrix prior TEI calculation.
    :param target_sum: Optional[float] (default: 1e6)
        After normalization, each observation (cell) has a total count equal to target_sum.
    :return:

    Example
    --------
    >>>
    """
    id_age_df_keep_subset, adata_counts, var_names_subset, sumx, sumx_recd, ps, psd =\
        get_psd(adata, gene_id, gene_age, keep, layer, normalize_total, log1p, target_sum)
    teimatrix = psd.dot(adata_counts.transpose()).transpose()
    pmatrix = sumx_recd.dot(teimatrix)
    tei = pmatrix.sum(1)
    phylostrata = list(set(id_age_df_keep_subset['Phylostrata']))
    pstrata_norm_by_sumx = np.empty((len(phylostrata), pmatrix.shape[0]))
    pstrata_norm_by_pmatrix_sum = np.empty((len(phylostrata), pmatrix.shape[0]))
    for pk_idx, pk in enumerate(phylostrata):
        pstrata_norm_by_sumx[pk_idx, ] = np.array(pmatrix[:, id_age_df_keep_subset['Phylostrata'].isin([pk])]
                                                  .sum(1)).flatten()
        pstrata_norm_by_pmatrix_sum[pk_idx, ] = np.array(pmatrix[:, id_age_df_keep_subset['Phylostrata']
                                                         .isin([pk])].sum(1) / tei).flatten()
    pstrata_norm_by_sumx_df = pd.DataFrame(pstrata_norm_by_sumx)
    pstrata_norm_by_sumx_df['ps'] = phylostrata
    pstrata_norm_by_sumx_df.set_index('ps', inplace=True)
    pstrata_norm_by_sumx_df.columns = adata.obs_names
    pstrata_norm_by_pmatrix_sum_df = pd.DataFrame(pstrata_norm_by_pmatrix_sum)
    pstrata_norm_by_pmatrix_sum_df['ps'] = phylostrata
    pstrata_norm_by_pmatrix_sum_df.set_index('ps', inplace=True)
    pstrata_norm_by_pmatrix_sum_df.columns = adata.obs_names
    if group_by is not None:
        pstrata_norm_by_sumx_df =\
            pstrata_norm_by_sumx_df.transpose().groupby(adata.obs[group_by]).mean().transpose()
        pstrata_norm_by_pmatrix_sum_df =\
            pstrata_norm_by_pmatrix_sum_df.transpose().groupby(adata.obs[group_by]).mean().transpose()
    if cumsum:
        pstrata_norm_by_sumx_df = pstrata_norm_by_sumx_df.cumsum(0)
        pstrata_norm_by_pmatrix_sum_df = pstrata_norm_by_pmatrix_sum_df.cumsum(0)
    if standard_scale is not None:
        if standard_scale == 0:
            pstrata_norm_by_sumx_df = pstrata_norm_by_sumx_df.apply(min_max_to_01, axis=1, raw=True)
            pstrata_norm_by_pmatrix_sum_df = pstrata_norm_by_pmatrix_sum_df.apply(min_max_to_01, axis=1, raw=True)
        if standard_scale == 1:
            pstrata_norm_by_sumx_df = pstrata_norm_by_sumx_df.apply(min_max_to_01, axis=0, raw=True)
            pstrata_norm_by_pmatrix_sum_df = pstrata_norm_by_pmatrix_sum_df.apply(min_max_to_01, axis=0, raw=True)
    return [pstrata_norm_by_sumx_df, pstrata_norm_by_pmatrix_sum_df]


def min_max_to_01(ndarray):
    """

    :param ndarray:
    :return:
    """
    ndarray_min = np.min(ndarray)
    ndarray_max = np.max(ndarray)
    if ndarray_min == ndarray_max:
        ndarray_min_max = [x - ndarray_min for x in ndarray]
    else:
        ndarray_min_max = [((x - ndarray_min)/(ndarray_max-ndarray_min)) for x in ndarray]
    return ndarray_min_max


def get_rematrix(adata, gene_id, gene_age, keep='min', layer=None, use=None, col_type='mean',
                 standard_scale=None, group_by=None, group_type='mean', normalize_total=False, log1p=False, target_sum=1e6):
    """
    This function computes relative expression profiles.

    In detail, if axis is set to None, the mean/median/sum expression is being computed over cells
    and, if group_by is not None, combined per given group by mean/median/sum.

    In detail, if axis is set to 0, the mean/median/sum relative expression profile is being computed over cells
    and, if group_by is not None, combined per given group by mean/median/sum as follows:

    f_c = (e_c - e_min)/(e_max - e_min)

    where e_min and e_max denote either the minimum/maximum mean/median/sum
    expression level over cells c.

    In detail, if axis is set to 1, the mean/median/sum relative expression profile is being computed over gene
    age classes (phylostrata) and, if group_by is not None, combined per given group by mean/median/sum as follows:

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

    :param adata: AnnData
        The annotated data matrix of shape n_obs × n_vars. Rows correspond to cells and columns to genes.
    :param gene_id: list
        Expects GeneID column from orthomap DataFrame.
    :param gene_age: list
        Expects Phylostratum column from orthomap DataFrame.
    :param keep: str (default: min)
        Either define 'min' (ascending pre-sorting) or 'max' (non-ascending pre-sorting) to keep duplicates.
    :param layer: Optional[str] (default: None)
        Layer to work on instead of X. If None, X is used.
    :param use: Optional[str] (default: None)
        Specify if counts from adata.X (default if None) should be combined per age group to calculate
        the relative expression or if the corresponding 'pmatrix' or 'teimatrix' should be used.
        If layer is not None adata.X refers to adata.layers[layer].
    :param col_type: str (default: mean)
        Specify how the counts should be combined, either by 'mean', 'median', 'sum', 'min' or 'max'.
    :param standard_scale: int (default: None)
        Whether or not to standardize the given axis (0: rows, gene age class, 1: columns, cells or group) between 0 and 1,
        meaning for each variable or group, subtract the minimum and divide each by its maximum.
    :param group_by:
    :param group_type: str (default: mean)
        Specify how the counts should be combined per group, either by 'mean', 'median', 'sum', 'min' or 'max'.
    :param normalize_total: bool (default: False)
        Normalize counts per cell prior TEI calculation.
    :param log1p: log1p: bool (default: False)
        Logarithmize the data matrix prior TEI calculation.
    :param target_sum: Optional[float] (default: 1e6)
        After normalization, each observation (cell) has a total count equal to target_sum.
    :return:

    Example
    --------
    >>>
    """
    id_age_df_keep_subset, adata_counts, var_names_subset, sumx, sumx_recd, ps, psd =\
        get_psd(adata, gene_id, gene_age, keep, layer, normalize_total, log1p, target_sum)
    teimatrix = psd.dot(adata_counts.transpose()).transpose()
    pmatrix = sumx_recd.dot(teimatrix)
    tei = pmatrix.sum(1)
    phylostrata = list(set(id_age_df_keep_subset['Phylostrata']))
    rematrix = np.empty((len(phylostrata), adata_counts.shape[0]))
    if use == 'pmatrix':
        for pk_idx, pk in enumerate(phylostrata):
            if col_type == 'mean':
                rematrix[pk_idx, ] = np.array(pmatrix[:, id_age_df_keep_subset['Phylostrata'].isin([pk])]
                                              .mean(1)).flatten()
            if col_type == 'median':
                rematrix[pk_idx, ] = np.apply_along_axis(
                    np.median, 1, pmatrix[:, id_age_df_keep_subset['Phylostrata'].isin([pk])].toarray()).flatten()
            if col_type == 'sum':
                rematrix[pk_idx, ] = np.array(pmatrix[:, id_age_df_keep_subset['Phylostrata'].isin([pk])]
                                              .sum(1)).flatten()
            if col_type == 'min':
                rematrix[pk_idx, ] = np.array(pmatrix[:, id_age_df_keep_subset['Phylostrata'].isin([pk])]
                                              .min(1).toarray()).flatten()
            if col_type == 'max':
                rematrix[pk_idx, ] = np.array(pmatrix[:, id_age_df_keep_subset['Phylostrata'].isin([pk])]
                                              .max(1).toarray()).flatten()
    if use == 'teimatrix':
        for pk_idx, pk in enumerate(phylostrata):
            if col_type == 'mean':
                rematrix[pk_idx, ] = np.array(teimatrix[:, id_age_df_keep_subset['Phylostrata'].isin([pk])]
                                              .mean(1)).flatten()
            if col_type == 'median':
                rematrix[pk_idx, ] = np.apply_along_axis(
                    np.median, 1, teimatrix[:, id_age_df_keep_subset['Phylostrata'].isin([pk])].toarray()).flatten()
            if col_type == 'sum':
                rematrix[pk_idx, ] = np.array(teimatrix[:, id_age_df_keep_subset['Phylostrata'].isin([pk])]
                                              .sum(1)).flatten()
            if col_type == 'min':
                rematrix[pk_idx, ] = np.array(teimatrix[:, id_age_df_keep_subset['Phylostrata'].isin([pk])]
                                              .min(1).toarray()).flatten()
            if col_type == 'max':
                rematrix[pk_idx, ] = np.array(teimatrix[:, id_age_df_keep_subset['Phylostrata'].isin([pk])]
                                              .max(1).toarray()).flatten()
    else:
        for pk_idx, pk in enumerate(phylostrata):
            if col_type == 'mean':
                rematrix[pk_idx, ] = np.array(adata_counts[:, id_age_df_keep_subset['Phylostrata'].isin([pk])]
                                              .mean(1)).flatten()
            if col_type == 'median':
                rematrix[pk_idx, ] = np.apply_along_axis(
                    np.median, 1, adata_counts[:, id_age_df_keep_subset['Phylostrata'].isin([pk])].toarray()).flatten()
            if col_type == 'sum':
                rematrix[pk_idx, ] = np.array(adata_counts[:, id_age_df_keep_subset['Phylostrata'].isin([pk])]
                                              .sum(1)).flatten()
            if col_type == 'min':
                rematrix[pk_idx, ] = np.array(adata_counts[:, id_age_df_keep_subset['Phylostrata'].isin([pk])]
                                              .min(1).toarray()).flatten()
            if col_type == 'max':
                rematrix[pk_idx, ] = np.array(adata_counts[:, id_age_df_keep_subset['Phylostrata'].isin([pk])]
                                              .max(1).toarray()).flatten()
    rematrix_df = pd.DataFrame(rematrix)
    rematrix_df['ps'] = phylostrata
    rematrix_df.set_index('ps', inplace=True)
    rematrix_df.columns = adata.obs_names
    if group_by is not None:
        if group_type == 'mean':
            rematrix_df = \
                rematrix_df.transpose().groupby(adata.obs[group_by]).mean().transpose()
        if group_type == 'median':
            rematrix_df = \
                rematrix_df.transpose().groupby(adata.obs[group_by]).median().transpose()
        if group_type == 'sum':
            rematrix_df = \
                rematrix_df.transpose().groupby(adata.obs[group_by]).sum().transpose()
        if group_type == 'min':
            rematrix_df = \
                rematrix_df.transpose().groupby(adata.obs[group_by]).min().transpose()
        if group_type == 'max':
            rematrix_df = \
                rematrix_df.transpose().groupby(adata.obs[group_by]).max().transpose()
    if standard_scale is not None:
        if standard_scale == 0:
            rematrix_df = rematrix_df.apply(min_max_to_01, axis=1, raw=True)
        if standard_scale == 1:
            rematrix_df = rematrix_df.apply(min_max_to_01, axis=0, raw=True)
    return rematrix_df
