#!/usr/bin/python
# -*- coding: UTF-8 -*-


"""
Author: Kristian K Ullrich
date: March 2023
email: ullrich@evolbio.mpg.de
License: GPL-3
"""


import os
import sys
import zipfile
import argparse
import pandas as pd
from orthomap import qlin
from ete3 import NCBITaxa


def define_parser():
    """
    A helper function for using `of2orthomap.py` via the terminal.

    :return: An argparse.ArgumentParser.

    :rtype: argparse.ArgumentParser
    """
    of2orthomap_example = '''of2orthomap example:

    # download OrthoFinder example:
    $ wget https://github.com/kullrich/orthomap/raw/main/examples/ensembl_105_orthofinder_Orthogroups.GeneCount.tsv.zip
    $ wget https://github.com/kullrich/orthomap/raw/main/examples/ensembl_105_orthofinder_Orthogroups.tsv.zip
    $ wget https://github.com/kullrich/orthomap/raw/main/examples/ensembl_105_orthofinder_species_list.tsv
    
    # extract orthomap:
    $ of2orthomap -seqname Danio_rerio.GRCz11.cds.longest -qt 7955 \\
      -sl ensembl_105_orthofinder_species_list.tsv \\
      -oc ensembl_105_orthofinder_Orthogroups.GeneCount.tsv.zip \\
      -og ensembl_105_orthofinder_Orthogroups.tsv.zip
    '''
    parser = argparse.ArgumentParser(
        prog='of2orthomap',
        usage='%(prog)s [options] [<arguments>...]',
        description='extract orthomap from orthofinder output for query species',
        epilog=of2orthomap_example,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    add_argparse_args(parser=parser)
    return parser


def add_argparse_args(parser: argparse.ArgumentParser):
    """
    This function attaches individual argument specifications to the parser.

    :param parser: An argparse.ArgumentParser.

    :type parser: argparse.ArgumentParser
    """
    parser.add_argument('-seqname', help='sequence name of the query species in OrthoFinder'
                                         '(see column names of  <Orthogroups.tsv>)')
    parser.add_argument('-qt', help='query species taxID (e.g. use <orthomap qlin -h> to get taxID)')
    parser.add_argument('-sl', help='species list as <OrthoFinder name><tab><species taxID> '
                                    '(only samples in this list will be processed)')
    parser.add_argument('-oc', help='specify OrthoFinder <Orthogroups.GeneCounts.tsv> (see Orthogroups directory)')
    parser.add_argument('-og', help='specify OrthoFinder <Orthogroups.tsv> (see Orthogroups directory)')
    parser.add_argument('-out', help='specify output file <orthomap.tsv> (default: orthomap.tsv)',
                        default='orthomap.tsv')
    parser.add_argument('-overwrite', help='specify if existing output file should be overwritten (default: True)',
                        default=True, type=bool)


def get_orthomap(seqname, qt, sl, oc, og, out=None, quiet=False, continuity=True, overwrite=True):
    """
    This function return an orthomap for a given query species and OrthoFinder input data.

    :param seqname: Sequence name of the query species used for OrthoFinder comparison.
    :param qt: Query species taxID.
    :param sl: Path to species list file containing <OrthoFinder name><tab><species taxID>.
    :param oc: Path to OrthoFinder result <Orthogroups.GeneCounts.tsv> file.
    :param og: Path to OrthoFinder result <Orthogroups.tsv> file.
    :param out: Path to output file.
    :param quiet: Specify if output should be quiet.
    :param continuity: Specify if continuity score should be calculated.
    :param overwrite: Specify if output should be overwritten.
    :return: A list of results such as:
    orthomap, species_list, younges_common_counts

    :type seqname: str
    :type qt: str
    :type sl: str
    :type oc: str
    :type og: str
    :type out: str
    :type quiet: bool
    :type continuity: bool
    :type overwrite: bool
    :rtype: list

    Example
    -------
    >>> from orthomap import datasets, of2orthomap, qlin
    >>> datasets.ensembl105(datapath='.')
    >>> query_orthomap = of2orthomap.get_orthomap(
    >>>     seqname='Danio_rerio.GRCz11.cds.longest',
    >>>     qt='7955',
    >>>     sl='ensembl_105_orthofinder_species_list.tsv',
    >>>     oc='ensembl_105_orthofinder_Orthogroups.GeneCount.tsv',
    >>>     og='ensembl_105_orthofinder_Orthogroups.tsv',
    >>>     out=None, quiet=False, continuity=True, overwrite=True)
    >>> query_orthomap
    """
    outhandle = None
    og_continuity_score = None
    ncbi = NCBITaxa()
    qname, qtid, qlineage, qlineagenames_dict, qlineagezip, qlineagenames, qlineagerev, qk = \
        qlin.get_qlin(qt=qt, quiet=True)
    query_lineage_topo = qlin.get_lineage_topo(qt)
    species_list = pd.read_csv(sl, sep='\t', header=None)
    species_list.columns = ['species', 'taxID']
    species_list['lineage'] = species_list.apply(lambda x: ncbi.get_lineage(x[1]), axis=1)
    species_list['youngest_common'] = [qlin.get_youngest_common(qlineage, x) for x in species_list.lineage]
    species_list['youngest_name'] = [list(x.values())[0] for x in [ncbi.get_taxid_translator([x])
                                                                   for x in list(species_list.youngest_common)]]
    if not quiet:
        print(seqname)
        print(qname)
        print(qt)
        print(species_list)
    youngest_common_counts_df = get_youngest_common_counts(qlineage, species_list)
    for node in query_lineage_topo.traverse('postorder'):
        nsplit = node.name.split('/')
        if len(nsplit) == 3:
            node.add_feature('species_count',
                             list(youngest_common_counts_df[youngest_common_counts_df.PStaxID.isin(
                                 [int(nsplit[1])])].counts)[0])
    oc_og_dict = {}
    continuity_dict = {}
    if os.path.basename(oc).split('.')[-1] == 'zip':
        oc_zip = zipfile.Path(oc, at='.'.join(os.path.basename(oc).split('.')[:-1]))
        oc_lines = oc_zip.open()
    else:
        oc_lines = open(oc, 'r')
    oc_species = next(oc_lines)
    if type(oc_species) == bytes:
        oc_species = oc_species.decode('utf-8').strip().split('\t')
    else:
        oc_species = oc_species.strip().split('\t')
    oc_qidx = [x for x, y in enumerate(oc_species) if y == seqname]
    if len(oc_qidx) == 0:
        print('\nError <-qname>: query species name not in orthofinder results, please check spelling\n'
              'e.g. <head -1 Orthogroups.GeneCounts.tsv>')
        sys.exit()
    for oc_line in oc_lines:
        if type(oc_line) == bytes:
            oc_og = oc_line.decode('utf-8').strip().split('\t')
        else:
            oc_og = oc_line.strip().split('\t')
        if int(oc_og[oc_qidx[0]]) == 0:
            continue
        if int(oc_og[oc_qidx[0]]) > 0:
            oc_og_hits = [oc_species[x+1] for x, y in enumerate(oc_og[1::][::-1][1::][::-1]) if int(y) > 0]
            # get list of the youngest common between query and all other species
            oc_og_hits_youngest_common = list(species_list.youngest_common[
                                                  [x for x, y in enumerate(species_list.species)
                                                   if y in oc_og_hits]])
            # evaluate all youngest common nodes to retain the oldest of them and assign as the orthogroup
            # ancestral state (gene age)
            if len(oc_og_hits_youngest_common) > 0:
                oc_og_oldest_common = qlin.get_oldest_common(qlineage, oc_og_hits_youngest_common)
                oc_og_dict[oc_og[0]] = oc_og_oldest_common
                if continuity:
                    continuity_dict[oc_og[0]] = get_youngest_common_counts(
                        qlineage, pd.DataFrame(oc_og_hits_youngest_common, columns=['youngest_common'])).counts
    oc_lines.close()
    if continuity:
        youngest_common_counts_df = youngest_common_counts_df.join(pd.DataFrame.from_dict(continuity_dict))
    omap = []
    if out:
        if os.path.exists(out) and not overwrite:
            print('\nError <-overwrite>: output file exists, please set to True if it should be overwritten\n')
            sys.exit()
        outhandle = open(out, 'w')
        if continuity:
            outhandle.write('seqID\tOrthogroup\tPSnum\tPStaxID\tPSname\tPScontinuity\n')
        else:
            outhandle.write('seqID\tOrthogroup\tPSnum\tPStaxID\tPSname\n')
    if os.path.basename(og).split('.')[-1] == 'zip':
        og_zip = zipfile.Path(og, at='.'.join(os.path.basename(og).split('.')[:-1]))
        og_lines = og_zip.open()
    else:
        og_lines = open(og, 'r')
    og_species = next(og_lines)
    if type(og_species) == bytes:
        og_species = og_species.decode('utf-8').strip().split('\t')
    else:
        og_species = og_species.strip().split('\t')
    og_qidx = [x for x, y in enumerate(og_species) if y == seqname]
    if len(oc_qidx) == 0:
        print('\nError <-qname>: query species name not in orthofinder results, please check spelling\n'
              'e.g. <head -1 Orthogroups.tsv>')
        sys.exit()
    for og_line in og_lines:
        if type(og_line) == bytes:
            og_og = og_line.decode('utf-8').strip().split('\t')
        else:
            og_og = og_line.strip().split('\t')
        if og_og[0] not in oc_og_dict:
            continue
        else:
            og_ps = qlineagenames[qlineagenames['PStaxID'] ==
                                  str(oc_og_dict[og_og[0]])].values.tolist()[0]
            og_ps_join = '\t'.join(og_ps)
            if continuity:
                og_continuity_score = get_continuity_score(og_og[0], youngest_common_counts_df)
            if out:
                if continuity:
                    [outhandle.write(x.replace(' ', '') + '\t' + og_og[0] + '\t' + og_ps_join + '\t' +
                                     str(og_continuity_score) + '\n') for x in og_og[og_qidx[0]].split(',')]
                else:
                    [outhandle.write(x.replace(' ', '') + '\t' + og_og[0] + '\t' + og_ps_join + '\n')
                     for x in og_og[og_qidx[0]].split(',')]
        if continuity:
            omap += [[x.replace(' ', ''), og_og[0], og_ps[0], og_ps[1], og_ps[2], og_continuity_score]
                     for x in og_og[og_qidx[0]].split(',')]
        else:
            omap += [[x.replace(' ', ''), og_og[0], og_ps[0], og_ps[1], og_ps[2]]
                     for x in og_og[og_qidx[0]].split(',')]
    og_lines.close()
    if out:
        outhandle.close()
    omap_df = pd.DataFrame(omap)
    if continuity:
        omap_df.columns = ['seqID', 'Orthogroup', 'PSnum', 'PStaxID', 'PSname', 'PScontinuity']
    else:
        omap_df.columns = ['seqID', 'Orthogroup', 'PSnum', 'PStaxID', 'PSname']
    omap_df['PSnum'] = [int(x) for x in list(omap_df['PSnum'])]
    return [omap_df, species_list, youngest_common_counts_df]


def get_counts_per_ps(omap_df, psnum_col='PSnum', pstaxid_col='PStaxID', psname_col='PSname'):
    """
    This function return counts per phylostratum.

    :param omap_df: DataFrame with orthomap results.
    :param psnum_col: Specify PSnum column name.
    :param pstaxid_col: Specify PStaxID column name.
    :param psname_col: Specify PSname column name.
    :return: DataFrame with counts per phylostratum.

    :type omap_df: pandas.DataFrame
    :type psnum_col: str
    :type pstaxid_col: str
    :type psname_col: str
    :rtype: pandas.DataFrame

    Example
    -------
    >>> from orthomap import datasets, of2orthomap, qlin
    >>> datasets.ensembl105(datapath='.')
    >>> query_orthomap = of2orthomap.get_orthomap(
    >>>     seqname='Danio_rerio.GRCz11.cds.longest',
    >>>     qt='7955',
    >>>     sl='ensembl_105_orthofinder_species_list.tsv',
    >>>     oc='ensembl_105_orthofinder_Orthogroups.GeneCount.tsv',
    >>>     og='ensembl_105_orthofinder_Orthogroups.tsv',
    >>>     out=None, quiet=False, continuity=True, overwrite=True)
    >>> of2orthomap.get_counts_per_ps(query_orthomap[0],
    >>>     psnum_col='PSnum',
    >>>     pstaxid_col='PStaxID',
    >>>     psname_col='PSname')
    """
    counts_df = pd.DataFrame(omap_df.value_counts(psnum_col))
    counts_df.columns = ['counts']
    counts_df.reset_index(inplace=True)
    if pstaxid_col:
        counts_df = counts_df.merge(omap_df[~omap_df.duplicated(psnum_col)][[psnum_col, pstaxid_col]], on=psnum_col)
    if psname_col:
        counts_df = counts_df.merge(omap_df[~omap_df.duplicated(psnum_col)][[psnum_col, psname_col]], on=psnum_col)
    counts_df.set_index(psnum_col, inplace=True, drop=False)
    counts_df = counts_df.sort_index()
    return counts_df


def get_youngest_common_counts(qlineage, species_list):
    """
    This function return LCA counts for a given query species lineage.

    :param qlineage: Query lineage information.
    :param species_list: Species list.
    :return: DataFrame with LCA counts.

    :type qlineage: list
    :type species_list: pandas.DataFrame
    :rtype: pandas.DataFrame

    Example
    -------
    >>> from orthomap import datasets, of2orthomap, qlin
    >>> datasets.ensembl105(datapath='.')
    """
    counts_df = pd.DataFrame(qlineage, columns=['lineage'])
    counts_df.set_index('lineage', inplace=True)
    counts_df = pd.concat([counts_df, species_list['youngest_common'].value_counts()], join='outer', axis=1)
    counts_df.columns = ['counts']
    counts_df['PStaxID'] = counts_df.index.values
    counts_df['PSnum'] = list(range(len(counts_df['PStaxID'])))
    return counts_df


def get_continuity_score(og_name, youngest_common_counts_df):
    """
    This function calculates a continuity score for a given orthologous group
    and its corresponding LCA counts.

    :param og_name: Orthologous group name.
    :param youngest_common_counts_df: DataFrame with LCA counts.
    :return: Continuity score.

    :type og_name: str
    :type youngest_common_counts_df: pandas.DataFrame
    :rtype: float

    Example
    -------
    >>>
    """
    og_continuity_score = 0.0
    og_df = youngest_common_counts_df[~youngest_common_counts_df['counts'].isna()][og_name]
    og_lca = (~og_df.isna()).idxmax()
    og_lca_pos = og_df.index.get_loc(og_lca)
    og_lca_df = og_df.iloc[og_lca_pos:]
    og_lca_df_counts = og_lca_df.isna().value_counts()
    if False in og_lca_df_counts:
        og_continuity_score = og_lca_df_counts[False] / len(og_lca_df)
    return og_continuity_score


def main():
    """
    The main function that is being called when `of2orthomap` is used via the terminal.
    """
    parser = define_parser()
    args = parser.parse_args()
    print(args)
    if not args.seqname:
        parser.print_help()
        print('\nError <-seqname>: Please specify query species name in OrthoFinder and taxID')
        sys.exit()
    if not args.qt:
        parser.print_help()
        print('\nError <-qt>: Please specify query species taxID')
        sys.exit()
    if not args.sl:
        parser.print_help()
        print('\nError <-sl>: Please specify species list as <OrthoFinder name><tab><species taxID>')
        sys.exit()
    if not args.oc:
        parser.print_help()
        print('\nError <-oc>: Please specify OrthoFinder <Orthogroups.GeneCounts.tsv> (see Orthogroups directory)')
        sys.exit()
    if not args.og:
        parser.print_help()
        print('\nError <-og>: Please specify OrthoFinder <Orthogroups.tsv> (see Orthogroups directory)')
        sys.exit()
    get_orthomap(seqname=args.seqname, qt=args.qt, sl=args.sl, oc=args.oc, og=args.og, out=args.out, quiet=False,
                 continuity=True, overwrite=args.overwrite)


if __name__ == '__main__':
    main()
