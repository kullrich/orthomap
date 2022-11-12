#!/usr/bin/python
# -*- coding: UTF-8 -*-


"""
Author: Kristian K Ullrich
date: October 2022
email: ullrich@evolbio.mpg.de
License: GPL-3
"""


import sys
import argparse
import pandas as pd
from orthomap import qlin
from ete3 import NCBITaxa


def define_parser():
    """

    :return:
    """
    of2orthomap_example = '''example:

    #
    of2orthomap -seqname -qt 10090 -sl -oc -og
    '''
    parser = argparse.ArgumentParser(prog='of2orthomap', usage='%(prog)s [options] [<arguments>...]',
                                     description='extract orthomap from orthofinder output for query species',
                                     epilog=of2orthomap_example,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    add_argparse_args(parser=parser)
    return parser


def add_argparse_args(parser: argparse.ArgumentParser):
    """

    :param parser:
    :return:
    """
    parser.add_argument('-seqname', help='sequence name of the query species in orthofinder'
                                         '(see column names of  <Orthogroups.tsv>)')
    parser.add_argument('-qt', help='query species taxid (e.g. use <orthomap qlin -h> to get taxid)')
    parser.add_argument('-sl', help='species list as <orthofinder name><tab><species taxid> '
                                    '(only samples in this list will be processed)')
    parser.add_argument('-oc', help='specify orthofinder <Orthogroups.GeneCounts.tsv> (see Orthogroups directory)')
    parser.add_argument('-og', help='specify orthofinder <Orthogroups.tsv> (see Orthogroups directory)')
    parser.add_argument('-out', help='specify output file <orthomap.tsv> (default: orthomap.tsv)',
                        default='orthomap.tsv')


def get_orthomap(seqname, qt, sl, oc, og, out=None, quite=False):
    """

    :param seqname:
    :param qt:
    :param sl:
    :param oc:
    :param og:
    :param out:
    :param quite:
    :return:
    """
    ncbi = NCBITaxa()
    qname, qtid, qlineage, qlineagenames_dict, qlineagezip, qlineagenames, qlineagerev, qk = qlin.get_qlin(qt=qt, quite=True)
    query_lineage_topo = qlin.get_lineage_topo(qt)
    species_list = pd.read_csv(sl, sep='\t', header=None)
    species_list.columns = ['species', 'taxID']
    species_list['lineage'] = species_list.apply(lambda x: ncbi.get_lineage(x[1]), axis=1)
    species_list['youngest_common'] = [qlin.get_youngest_common(qlineage, x) for x in species_list.lineage]
    species_list['youngest_name'] = [list(x.values())[0] for x in [ncbi.get_taxid_translator([x])
                                                                   for x in list(species_list.youngest_common)]]
    if not quite:
        print(seqname)
        print(qname)
        print(qt)
        print(species_list)
    oc_og_dict = {}
    with open(oc, 'r') as oc_lines:
        oc_species = next(oc_lines).strip().split('\t')
        oc_qidx = [x for x, y in enumerate(oc_species) if y == seqname]
        if len(oc_qidx) == 0:
            print('\nError <-qname>: query species name not in orthofinder results, please check spelling\n'
                  'e.g. <head -1 Orthogroups.GeneCounts.tsv>')
            sys.exit()
        for oc_line in oc_lines:
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
    omap = []
    if out:
        outhandle = open(out, 'w')
        outhandle.write('seqID\tOrthogroup\tPSnum\tPStaxID\tPSname\n')
    with open(og, 'r') as og_lines:
        og_species = next(og_lines).strip().split('\t')
        og_qidx = [x for x, y in enumerate(og_species) if y == seqname]
        if len(oc_qidx) == 0:
            print('\nError <-qname>: query species name not in orthofinder results, please check spelling\n'
                  'e.g. <head -1 Orthogroups.tsv>')
            sys.exit()
        for og_line in og_lines:
            og_og = og_line.strip().split('\t')
            if og_og[0] not in oc_og_dict:
                continue
            else:
                og_ps = qlineagenames[qlineagenames['PStaxID'] ==
                                      str(oc_og_dict[og_og[0]])].values.tolist()[0]
                og_ps_join = '\t'.join(og_ps)
                if out:
                    [outhandle.write(x.replace(' ', '') + '\t' + og_og[0] + '\t' + og_ps_join + '\n')
                     for x in og_og[og_qidx[0]].split(',')]
            omap += [[x.replace(' ', ''), og_og[0], og_ps[0], og_ps[1], og_ps[2]]
                     for x in og_og[og_qidx[0]].split(',')]
    if out:
        outhandle.close()
    omap_df = pd.DataFrame(omap)
    omap_df.columns = ['seqID', 'Orthogroup', 'PSnum', 'PStaxID', 'PSname']
    omap_df['PSnum'] = [int(x) for x in list(omap_df['PSnum'])]
    return [omap_df, species_list]


def get_counts_per_ps(omap_df, psnum_col='PSnum', pstaxid_col='PStaxID', psname_col='PSname'):
    counts_df = pd.DataFrame(omap_df[psnum_col].value_counts())
    counts_df.columns = ['counts']
    counts_df[psnum_col] = list(list(omap_df[psnum_col].value_counts().index.values))
    if pstaxid_col:
        counts_df[pstaxid_col] = list(list(omap_df[pstaxid_col].value_counts().index.values))
    if psname_col:
        counts_df[psname_col] = list(list(omap_df[psname_col].value_counts().index.values))
    counts_df = counts_df.sort_index()
    return counts_df


def main():
    """

    :return:
    """
    parser = define_parser()
    args = parser.parse_args()
    print(args)
    if not args.seqname:
        parser.print_help()
        print('\nError <-seqname>: Please specify query species name in orthofinder and taxid')
        sys.exit()
    if not args.qt:
        parser.print_help()
        print('\nError <-qt>: Please specify query species taxid')
        sys.exit()
    if not args.sl:
        parser.print_help()
        print('\nError <-sl>: Please specify species list as <orthofinder name><tab><species taxid>')
        sys.exit()
    if not args.oc:
        parser.print_help()
        print('\nError <-oc>: Please specify orthofinder <Orthogroups.GeneCounts.tsv> (see Orthogroups directory)')
        sys.exit()
    if not args.og:
        parser.print_help()
        print('\nError <-og>: Please specify orthofinder <Orthogroups.tsv> (see Orthogroups directory)')
        sys.exit()
    get_orthomap(args.qname, args.qt, args.sl, args.oc, args.og, args.out)


if __name__ == '__main__':
    main()