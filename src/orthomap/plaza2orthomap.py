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
import argparse
import pandas as pd
import numpy as np
from orthomap import of2orthomap, qlin
from ete3 import NCBITaxa


def define_parser():
    """
    A helper function for using `plaza2orthomap.py` via the terminal.

    :return: An argparse.ArgumentParser.

    :rtype: argparse.ArgumentParser
    """
    plaza2orthomap_example = '''example:

    # using Orthologous gene family 
    $ plaza2orthomap -qt 3702 -og genefamily_data.ORTHOFAM.csv -sl species_information.csv -out 3702.orthofam.orthomap
    
    # using Homologous gene family 
    $ plaza2orthomap -qt 3702 -og genefamily_data.HOMFAM.csv -sl species_information.csv -out 3702.homfam.orthomap
    '''
    parser = argparse.ArgumentParser(
        prog='plaza2orthomap',
        usage='%(prog)s [options] [<arguments>...]',
        description='extract orthomap from plaza gene family data for query species',
        epilog=plaza2orthomap_example,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    add_argparse_args(parser=parser)
    return parser


def add_argparse_args(parser: argparse.ArgumentParser):
    """
    This function attaches individual argument specifications to the parser.

    :param parser: An argparse.ArgumentParser.

    :type parser: argparse.ArgumentParser
    """
    parser.add_argument('-qt',
                        help='query species taxID (e.g. use <orthomap qlin -h> to get taxID)')
    parser.add_argument('-og',
                        help='specify plaza gene family file <genefamily_data.ORTHOFAM.csv> or '
                             'genefamily_data.HOMFAM.csv')
    parser.add_argument('-sl',
                        help='specify plaza species information file <species_information.csv>')
    parser.add_argument('-out',
                        help='specify output file <orthomap.tsv> (default: orthomap.tsv)',
                        default='orthomap.tsv')
    parser.add_argument('-overwrite',
                        help='specify if existing output file should be overwritten (default: True)',
                        default=True,
                        type=bool)


def _get_species_tax_id(species_name_list, species_list):
    """
    A helper function to map plaza species short name and species taxID.

    :param species_name_list: List of plaza species short names.
    :param species_list: DataFrame with plaza species information.
    :return: List of species taxID.

    :type species_name_list: list
    :type species_list: pandas.DataFrame
    :rtype: list
    """
    species_taxid_list = []
    for species_name in species_name_list:
        species_taxid = list(species_list['tax_id'][species_list['species'] == species_name])
        if len(species_taxid) > 0:
            species_taxid_list.append(species_taxid[0])
    return species_taxid_list


def get_plaza_orthomap(qt,
                       og,
                       sl,
                       out=None,
                       quiet=False,
                       continuity=True,
                       overwrite=True):
    """
    This function return an orthomap for a given query species and plaza gene family data.

    :param qt: Query species taxID.
    :param og: Path to plaza gene families <genefamily_data.ORTHOFAM.csv> file.
    :param sl: Path to plaza species information <species_information.csv> file.
    :param out: Path to output file.
    :param quiet: Specify if output should be quiet.
    :param continuity: Specify if continuity score should be calculated.
    :param overwrite: Specify if output should be overwritten.
    :return: A list of results such as:
             orthomap, species_list, youngest_common_counts

    :type qt: str
    :type og: str
    :type sl: str
    :type out: str
    :type quiet: bool
    :type continuity: bool
    :type overwrite: bool
    :rtype: list

    Example
    -------
    >>>
    """
    outhandle = None
    og_continuity_score = None
    ncbi = NCBITaxa()
    qname,\
        qtid,\
        qlineage,\
        qlineagenames_dict,\
        qlineagezip,\
        qlineagenames,\
        qlineagerev,\
        qk = qlin.get_qlin(qt=qt,
                           quiet=True)
    query_lineage_topo = qlin.get_lineage_topo(qt)
    species_list = pd.read_csv(sl, sep='\t', header=None, comment='#')
    species_list.columns = ['species', 'common_name', 'tax_id', 'source', 'data_provider', 'pubmed_id']
    qt_species = list(species_list['species'][species_list['tax_id'] == int(qt)])
    if len(qt_species) == 0:
        print('\nError <-qt>: query species taxID not in plaza results, please check taxID.')
        sys.exit()
    ogs = pd.DataFrame(pd.read_csv(og, sep='\t', header=None, comment='#'))
    ogs.columns = ['gf_id', 'species', 'gene_id']
    ogs_grouped = ogs.groupby('gf_id')['species'].apply(set).apply(list).apply(_get_species_tax_id,
                                                                               species_list=species_list)
    ogs_grouped_qt = pd.DataFrame(ogs_grouped[ogs_grouped.apply(lambda x: int(qt) in x)])
    ogs_qt = ogs[ogs['gf_id'].isin(ogs_grouped_qt.index)]
    ogs_qt_red = ogs_qt[ogs_qt['species'].isin(qt_species)]
    ogs_qt_red_grouped = ogs_qt_red.groupby('gf_id')['gene_id'].apply(list)
    ogs_grouped_qt['gene_id'] = ogs_qt_red_grouped
    ogs_grouped_qt_species = np.sort(list(set([x[0] for x in ogs_grouped_qt['species'].to_dict().values()])))
    ogs_grouped_qt_species_names = [qlin.get_qlin(qt=x,
                                                  quiet=True)[0] for x in ogs_grouped_qt_species]
    species_list_df = pd.DataFrame(ogs_grouped_qt_species_names,
                                   columns=['species'])
    species_list_df['taxID'] = ogs_grouped_qt_species
    species_list_df['lineage'] = species_list_df.apply(lambda x: ncbi.get_lineage(x[1]),
                                                       axis=1)
    species_list_df['youngest_common'] = [qlin.get_youngest_common(qlineage,
                                                                   x) for x in species_list_df.lineage]
    species_list_df['youngest_name'] = [list(x.values())[0] for x in [ncbi.get_taxid_translator([x])
                                                                      for x in list(species_list_df.youngest_common)]]
    if not quiet:
        print(qname)
        print(qt)
        print(species_list_df)
    youngest_common_counts_df = of2orthomap.get_youngest_common_counts(qlineage,
                                                                       species_list_df)
    for node in query_lineage_topo.traverse('postorder'):
        nsplit = node.name.split('/')
        if len(nsplit) == 3:
            node.add_feature('species_count',
                             list(youngest_common_counts_df[youngest_common_counts_df.PStaxID.isin(
                                 [int(nsplit[1])])].counts)[0])
    og_dict = {}
    continuity_dict = {}
    for og in ogs_grouped_qt.index:
        og_hits = np.sort(
            list(set([x[0] for x in ogs_grouped_qt[ogs_grouped_qt.index.isin([og])]['species'].to_dict().values()])))
        # get list of the youngest common between query and all other species
        og_hits_youngest_common = list(species_list_df.youngest_common[
                                           [x for x, y in enumerate(species_list_df.taxID)
                                            if y in og_hits]])
        # evaluate all youngest common nodes to retain the oldest of them and assign as the orthogroup
        # ancestral state (gene age)
        if len(og_hits_youngest_common) > 0:
            og_oldest_common = qlin.get_oldest_common(qlineage,
                                                      og_hits_youngest_common)
            og_dict[og] = og_oldest_common
            if continuity:
                continuity_dict[og] = \
                    of2orthomap.get_youngest_common_counts(qlineage,
                                                           pd.DataFrame(og_hits_youngest_common,
                                                                        columns=['youngest_common'])).counts
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
    for og in ogs_grouped_qt.index:
        og_tmp = ogs_grouped_qt[ogs_grouped_qt.index.isin([og])]
        og_ps = qlineagenames[qlineagenames['PStaxID'] ==
                              str(og_dict[og])].values.tolist()[0]
        og_ps_join = '\t'.join(og_ps)
        if continuity:
            og_continuity_score = of2orthomap.get_continuity_score(og,
                                                                   youngest_common_counts_df)
            if out:
                if continuity:
                    [outhandle.write(x.replace(' ', '') + '\t' + og + '\t' + og_ps_join + '\t' +
                                     str(og_continuity_score) + '\n') for x in list(og_tmp['gene_id'])[0]]
                else:
                    [outhandle.write(x.replace(' ', '') + '\t' + og + '\t' + og_ps_join + '\n')
                     for x in list(og_tmp['gene_id'])[0]]
        if continuity:
            omap += [[x.replace(' ', ''), og, og_ps[0], og_ps[1], og_ps[2], og_continuity_score]
                     for x in list(og_tmp['gene_id'])[0]]
        else:
            omap += [[x.replace(' ', ''), og, og_ps[0], og_ps[1], og_ps[2]]
                     for x in list(og_tmp['gene_id'])[0]]
    if out:
        outhandle.close()
    omap_df = pd.DataFrame(omap)
    if continuity:
        omap_df.columns = ['seqID',
                           'Orthogroup',
                           'PSnum',
                           'PStaxID',
                           'PSname',
                           'PScontinuity']
    else:
        omap_df.columns = ['seqID',
                           'Orthogroup',
                           'PSnum',
                           'PStaxID',
                           'PSname']
    omap_df['PSnum'] = [int(x) for x in list(omap_df['PSnum'])]
    return [omap_df,
            species_list_df,
            youngest_common_counts_df]


def main():
    """
    The main function that is being called when `plaza2orthomap` is used via the terminal.
    """
    parser = define_parser()
    args = parser.parse_args()
    print(args)
    if not args.qt:
        parser.print_help()
        print('\nError <-qt>: Please specify query species taxID')
        sys.exit()
    if not args.og:
        parser.print_help()
        print('\nError <-og>: Please specify plaza gene family file <genefamily_data.ORTHOFAM.csv> or '
              '<genefamily_data.HOMFAM.csv>')
        sys.exit()
    if not args.sl:
        parser.print_help()
        print('\nError <-sl>: Please specify plaza species information file <species_information.csv>')
        sys.exit()
    get_plaza_orthomap(qt=args.qt,
                       og=args.og,
                       sl=args.sl,
                       out=args.out,
                       overwrite=args.overwrite)


if __name__ == '__main__':
    main()
