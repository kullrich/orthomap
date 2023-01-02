#!/usr/bin/python
# -*- coding: UTF-8 -*-


"""
Author: Kristian K Ullrich
date: November 2022
email: ullrich@evolbio.mpg.de
License: GPL-3
"""


import argparse
import os
import sys

import pandas as pd
from ete3 import NCBITaxa

from orthomap import of2orthomap, qlin


def define_parser():
    """
    A helper function for using `eggnog2orthomap.py` via the terminal.

    :return: argparse.ArgumentParser
    """
    eggnog2orthomap_example = '''example:

    #
    $ eggnog2orthomap -qt 10090 -og
    '''
    parser = argparse.ArgumentParser(
        prog='eggnog2orthomap',
        usage='%(prog)s [options] [<arguments>...]',
        description='extract orthomap from eggnog output for query species',
        epilog=eggnog2orthomap_example,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    add_argparse_args(parser=parser)
    return parser


def add_argparse_args(parser: argparse.ArgumentParser):
    """
    This function attaches individual argument specifications to the parser.

    :param parser: argparse.ArgumentParser
    """
    parser.add_argument('-qt', help='query species taxid (e.g. use <orthomap qlin -h> to get taxid)')
    parser.add_argument('-og', help='specify eggnog <e6.og2seqs_and_species.tsv>')
    parser.add_argument('-subset', help='specify file of orthologous groups to include'
                                        '<e6.og2parents_and_children.new.tsv>')
    parser.add_argument('-out', help='specify output file <orthomap.tsv> (default: orthomap.tsv)',
                        default='orthomap.tsv')
    parser.add_argument('-overwrite', help='specify if existing output file should be overwritten (default: True)',
                        default=True, type=bool)


def get_eggnog_orthomap(qt, og, subset=None, out=None, quite=False, continuity=True, overwrite=True):
    """

    :param qt:
    :param og:
    :param subset:
    :param out:
    :param quite:
    :param continuity:
    :return:
    """
    ncbi = NCBITaxa()
    qname, qtid, qlineage, qlineagenames_dict, qlineagezip, qlineagenames, qlineagerev, qk = qlin.get_qlin(qt=qt, quite=True)
    query_lineage_topo = qlin.get_lineage_topo(qt)
    if subset is not None:
        subset_dict = {}
        with open(subset, 'r') as subset_ogs:
            for subset_tmp in subset_ogs:
                sog_name = subset_tmp.strip().split('\t')[0]
                subset_dict[sog_name] = []
    ogs_dict = {}
    species_list = []
    with open(og, 'r') as ogs:
        for og_line in ogs:
            col1_taxonomic_level, col2_og_name, col3_number_of_species, col4_number_of_members,\
            col5_comma_separated_list_of_species, col6_comma_separated_list_of_members = og_line.strip().split('\t')
            col5_comma_separated_list_of_species = col5_comma_separated_list_of_species.split(',')
            if subset is not None:
                if col2_og_name not in subset_dict:
                    continue
                else:
                    if str(qtid) in col5_comma_separated_list_of_species:
                        col6_comma_separated_list_of_members = col6_comma_separated_list_of_members.split(',')
                        q_genes = [x for x in col6_comma_separated_list_of_members if x.split('.')[0] == str(qtid)]
                        ogs_dict[col2_og_name] = [col2_og_name, col5_comma_separated_list_of_species, q_genes]
                        species_list += col5_comma_separated_list_of_species
            else:
                if str(qtid) in col5_comma_separated_list_of_species:
                    col6_comma_separated_list_of_members = col6_comma_separated_list_of_members.split(',')
                    q_genes = [x for x in col6_comma_separated_list_of_members if x.split('.')[0] == str(qtid)]
                    ogs_dict[col2_og_name] = [col2_og_name, col5_comma_separated_list_of_species, q_genes]
                    species_list += col5_comma_separated_list_of_species
    species_list = list(set(species_list))
    if len(species_list) == 0:
        print('\nError <-qt>: query species taxID not in eggnog results, please check taxID.')
        sys.exit()
    species_names = [qlin.get_qlin(qt=x, quite=True)[0] for x in species_list]
    species_list_df = pd.DataFrame(species_names, columns=['species'])
    species_list_df['taxID'] = [int(x) for x in species_list]
    species_list_df['lineage'] = species_list_df.apply(lambda x: ncbi.get_lineage(x[1]), axis=1)
    species_list_df['youngest_common'] = [qlin.get_youngest_common(qlineage, x) for x in species_list_df.lineage]
    species_list_df['youngest_name'] = [list(x.values())[0] for x in [ncbi.get_taxid_translator([x])
                                                                      for x in list(species_list_df.youngest_common)]]
    if not quite:
        print(qname)
        print(qt)
        print(species_list_df)
    youngest_common_counts_df = of2orthomap.get_youngest_common_counts(qlineage, species_list_df)
    for node in query_lineage_topo.traverse('postorder'):
        nsplit = node.name.split('/')
        if len(nsplit) == 3:
            node.add_feature('species_count',
                             list(youngest_common_counts_df[youngest_common_counts_df.PStaxID.isin(
                                 [int(nsplit[1])])].counts)[0])
    og_dict = {}
    continuity_dict = {}
    for og in ogs_dict.keys():
        og_tmp = ogs_dict[og]
        og_hits = [int(x) for x in og_tmp[1]]
        # get list of the youngest common between query and all other species
        og_hits_youngest_common = list(species_list_df.youngest_common[
                                            [x for x, y in enumerate(species_list_df.taxID)
                                            if y in og_hits]])
        # evaluate all youngest common nodes to retain the oldest of them and assign as the orthogroup
        # ancestral state (gene age)
        if len(og_hits_youngest_common) > 0:
            og_oldest_common = qlin.get_oldest_common(qlineage, og_hits_youngest_common)
            og_dict[og_tmp[0]] = og_oldest_common
            if continuity:
                continuity_dict[og_tmp[0]] =\
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
    for og in ogs_dict.keys():
        og_tmp = ogs_dict[og]
        og_ps = qlineagenames[qlineagenames['PStaxID'] ==
                                str(og_dict[og_tmp[0]])].values.tolist()[0]
        og_ps_join = '\t'.join(og_ps)
        if continuity:
            og_continuity_score = of2orthomap.get_continuity_score(og_tmp[0], youngest_common_counts_df)
            if out:
                if continuity:
                    [outhandle.write(x.replace(' ', '') + '\t' + og_tmp[0] + '\t' + og_ps_join + '\t' +
                                     str(og_continuity_score) + '\n') for x in og_tmp[2]]
                else:
                    [outhandle.write(x.replace(' ', '') + '\t' + og_tmp[0] + '\t' + og_ps_join + '\n')
                     for x in og_tmp[2]]
        if continuity:
            omap += [[x.replace(' ', ''), og_tmp[0], og_ps[0], og_ps[1], og_ps[2], og_continuity_score]
                     for x in og_tmp[2]]
        else:
            omap += [[x.replace(' ', ''), og_tmp[0], og_ps[0], og_ps[1], og_ps[2]]
                     for x in og_tmp[2]]
    if out:
        outhandle.close()
    omap_df = pd.DataFrame(omap)
    if continuity:
        omap_df.columns = ['seqID', 'Orthogroup', 'PSnum', 'PStaxID', 'PSname', 'PScontinuity']
    else:
        omap_df.columns = ['seqID', 'Orthogroup', 'PSnum', 'PStaxID', 'PSname']
    omap_df['PSnum'] = [int(x) for x in list(omap_df['PSnum'])]
    return [omap_df, species_list_df, youngest_common_counts_df]


def main():
    """
    The main function that is being called when `eggnog2orthomap.py` is used via the terminal.
    """
    parser = define_parser()
    args = parser.parse_args()
    print(args)
    if not args.qt:
        parser.print_help()
        print('\nError <-qt>: Please specify query species taxid')
        sys.exit()
    if not args.og:
        parser.print_help()
        print('\nError <-og>: Please specify eggnog <e6.og2seqs_and_species.tsv>')
        sys.exit()
    get_eggnog_orthomap(args.qt, args.og, subset=args.subset, out=args.out, overwrite=args.overwrite)


if __name__ == '__main__':
    main()
