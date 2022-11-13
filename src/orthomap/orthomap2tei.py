"""
Author: Kristian K Ullrich
date: October 2022
email: ullrich@evolbio.mpg.de
License: GPL-3
"""


import scipy
import numpy as np
import pandas as pd
import anndata as ad
from alive_progress import alive_bar


def read_orthomap(orthomapfile):
    """

    :param orthomapfile:
    :return:
    """
    orthomap = pd.read_csv(orthomapfile, delimiter='\t')
    return orthomap


def geneset_overlap(geneset1, geneset2):
    """

    :param geneset1:
    :param geneset2:
    :return:
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
    """
    replace_dict = {}
    for i, j in enumerate(xmatch):
        replace_dict[j] = xreplace[i]
    x_new = [replace_dict[x] for x in x_orig]
    return x_new


def keep_min_max(df, keep='min', dup_col=['GeneID'], sort_col=['Phylostrata']):
    """

    :param df:
    :param keep:
    :param dup_col:
    :param sort_col:
    :return:
    """
    if keep == 'min':
        df_sorted = df.sort_values(by=sort_col, ascending=True)
    if keep == 'max':
        df_sorted = df.sort_values(by=sort_col, ascending=False)
    keep_idx = ~df_sorted.duplicated(dup_col)
    df_out = df_sorted[keep_idx]
    return df_out


def get_tei(adata, gene_id, gene_age, keep='min', add=True, boot=False, bt=10):
    """
    This function computes the phylogenetically based transcriptome evolutionary
    index (TEI) similar to Domazet-Loso & Tautz, 2010.

    If the parameter boot is set to true,
    the strata values are sampled and the global TEI
    is calculated bt times.

    :param adata:
    :param gene_id:
    :param gene_age:
    :param keep:
    :param add:
    :param boot:
    :param bt:
    :return:
    """
    id_age_df = pd.DataFrame(data={'GeneID': gene_id, 'Phylostrata': gene_age})
    # check and drop duplicated GeneID
    id_age_df_keep = keep_min_max(id_age_df, keep=keep, dup_col=['GeneID'], sort_col=['Phylostrata'])
    # get overlap
    gene_intersection = pd.Index(id_age_df_keep['GeneID']).intersection(adata.var_names)
    id_age_df_keep_subset = id_age_df_keep.loc[id_age_df_keep['GeneID'].isin(gene_intersection)]
    id_age_df_keep_subset = id_age_df_keep_subset.sort_values('GeneID')
    adata_counts = adata.X
    adata_counts = adata_counts[:, adata.var_names.isin(id_age_df_keep_subset['GeneID'])]
    var_names_subset = adata.var_names[adata.var_names.isin(id_age_df_keep_subset['GeneID'])]
    var_names_subset_idx = var_names_subset.sort_values(return_indexer=True)[1]
    adata_counts = adata_counts[:, var_names_subset_idx]
    sumx = adata_counts.sum(1)
    ps = np.array(id_age_df_keep_subset['Phylostrata'])
    psd = scipy.sparse.diags(ps)
    teisum = psd.dot(adata_counts.transpose()).transpose().sum(1)
    tei = teisum/sumx
    tei_df = pd.DataFrame(tei, columns=['tei'])
    tei_df.index=adata.obs_names
    if add:
        adata.obs['tei'] = tei_df
    if boot:
        with alive_bar(bt) as bar:
            for i in range(bt):
                np.random.shuffle(ps)
                psd = scipy.sparse.diags(ps)
                teisum = psd.dot(adata_counts.transpose()).transpose().sum(1)
                tei = teisum / sumx
                tei_df[i] = tei
                bar()
    return tei_df


def get_pmatrix(adata, gene_id, gene_age, keep='min', add_obs=True, add_var=True):
    """
    This function computes the partial transcriptome evolutionary index (TEI) values for each single gene.

    In detail, each gene gets a TEI contribution profile as follows:
    \deqn{TEI_is = f_is * ps_i}
    where TEI_is is the partial TEI value of gene i,
    \eqn{f_is = e_is / \sum e_is} and \eqn{ps_i} is the phylostratum of gene i.

    The partial TEI matrix can be used to perform different cluster
    analyses and also gives an overall impression of the contribution of each
    gene to the global TEI pattern.

    :param adata:
    :param gene_id:
    :param gene_age:
    :param keep:
    :param add_obs:
    :param add_var:
    :return:
    """
    id_age_df = pd.DataFrame(data={'GeneID': gene_id, 'Phylostrata': gene_age})
    # check and drop duplicated GeneID
    id_age_df_keep = keep_min_max(id_age_df, keep=keep, dup_col=['GeneID'], sort_col=['Phylostrata'])
    # get overlap
    gene_intersection = pd.Index(id_age_df_keep['GeneID']).intersection(adata.var_names)
    id_age_df_keep_subset = id_age_df_keep.loc[id_age_df_keep['GeneID'].isin(gene_intersection)]
    id_age_df_keep_subset = id_age_df_keep_subset.sort_values('GeneID')
    adata_counts = adata.X
    adata_counts = adata_counts[:, adata.var_names.isin(id_age_df_keep_subset['GeneID'])]
    var_names_subset = adata.var_names[adata.var_names.isin(id_age_df_keep_subset['GeneID'])]
    var_names_subset_idx = var_names_subset.sort_values(return_indexer=True)[1]
    adata_counts = adata_counts[:, var_names_subset_idx]
    ps = np.array(id_age_df_keep_subset['Phylostrata'])
    psd = scipy.sparse.diags(ps)
    pmatrix = psd.dot(adata_counts.transpose()).transpose()
    adata_pmatrix = ad.AnnData(pmatrix)
    adata_pmatrix.obs_names = adata.obs_names
    adata_pmatrix.var_names = var_names_subset
    if add_var:
        for kv in adata.var.keys():
            adata_pmatrix.var[kv] = adata.var[kv][adata.var_names.isin(id_age_df_keep_subset['GeneID'])]
    if add_obs:
        for ko in adata.obs.keys():
            adata_pmatrix.obs[ko] = adata.obs[ko]
    return adata_pmatrix

def get_pstrata(adata, gene_id, gene_age, keep='min', cumsum=False, group_by=None):
    """
    This function computes the partial transcriptome evolutionary index (TEI) values combined for each strata.

    In detail, each gene gets a TEI contribution profile as follows:
    \deqn{TEI_is = f_is * ps_i}
    where TEI_is is the partial TEI value of gene i,
    \eqn{f_is = e_is / \sum e_is} and \eqn{ps_i} is the phylostratum of gene i.

    \eqn{TEI_is} values are combined per \eqn{ps}.

    The partial TEI values combined per strata give an overall impression
    of the contribution of each strata to the global TEI pattern.

    :param adata:
    :param gene_id:
    :param gene_age:
    :param keep:
    :param cumsum:
    :param group_by:
    :return:
    """
    id_age_df = pd.DataFrame(data={'GeneID': gene_id, 'Phylostrata': gene_age})
    # check and drop duplicated GeneID
    id_age_df_keep = keep_min_max(id_age_df, keep=keep, dup_col=['GeneID'], sort_col=['Phylostrata'])
    # get overlap
    gene_intersection = pd.Index(id_age_df_keep['GeneID']).intersection(adata.var_names)
    id_age_df_keep_subset = id_age_df_keep.loc[id_age_df_keep['GeneID'].isin(gene_intersection)]
    id_age_df_keep_subset = id_age_df_keep_subset.sort_values('GeneID')
    adata_counts = adata.X
    adata_counts = adata_counts[:, adata.var_names.isin(id_age_df_keep_subset['GeneID'])]
    var_names_subset = adata.var_names[adata.var_names.isin(id_age_df_keep_subset['GeneID'])]
    var_names_subset_idx = var_names_subset.sort_values(return_indexer=True)[1]
    adata_counts = adata_counts[:, var_names_subset_idx]
    sumx = adata_counts.sum(1)
    ps = np.array(id_age_df_keep_subset['Phylostrata'])
    psd = scipy.sparse.diags(ps)
    pmatrix = psd.dot(adata_counts.transpose()).transpose()
    pmatrix_sum = pmatrix.sum(1)
    phylostrata = list(set(id_age_df_keep_subset['Phylostrata']))
    pstrata_norm_by_sumx = np.empty((0, pmatrix.shape[0]))
    pstrata_norm_by_pmatrix_sum = np.empty((0, pmatrix.shape[0]))
    for pk in phylostrata:
        pstrata_norm_by_sumx = np.concatenate((pstrata_norm_by_sumx,
                                               (pmatrix[:, id_age_df_keep_subset['Phylostrata'].isin([pk])].sum(1) / sumx).flatten()))
        pstrata_norm_by_pmatrix_sum = np.concatenate((pstrata_norm_by_pmatrix_sum,
                                               (pmatrix[:, id_age_df_keep_subset['Phylostrata'].isin([pk])].sum(1) / pmatrix_sum).flatten()))
    pstrata_norm_by_sumx_df = pd.DataFrame(pstrata_norm_by_sumx)
    pstrata_norm_by_sumx_df['ps'] = phylostrata
    pstrata_norm_by_sumx_df.set_index('ps', inplace=True)
    pstrata_norm_by_sumx_df.columns = adata.obs_names
    pstrata_norm_by_pmatrix_sum_df = pd.DataFrame(pstrata_norm_by_pmatrix_sum)
    pstrata_norm_by_pmatrix_sum_df['ps'] = phylostrata
    pstrata_norm_by_pmatrix_sum_df.set_index('ps', inplace=True)
    pstrata_norm_by_pmatrix_sum_df.columns = adata.obs_names
    if group_by is not None:
        return
    if cumsum:
        pstrata_norm_by_sumx_df = pstrata_norm_by_sumx_df.cumsum(0)
        pstrata_norm_by_pmatrix_sum_df = pstrata_norm_by_pmatrix_sum_df.cumsum(0)
    return [pstrata_norm_by_sumx, pstrata_norm_by_pmatrix_sum]

def get_rematrix():
    return