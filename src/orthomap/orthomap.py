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
from ete3 import NCBITaxa


def define_parser():
    """

    :return:
    """
    orthomap_example = '''example:

    #
    orthomap -qname -qt 10090 -sl -oc -og
    '''
    parser = argparse.ArgumentParser(prog='orthomap', usage='%(prog)s [options] [<arguments>...]',
                                     description='extract orthomap from orthofinder output for query species',
                                     epilog=orthomap_example,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    add_argparse_args(parser=parser)
    return parser


def add_argparse_args(parser: argparse.ArgumentParser):
    """

    :param parser:
    :return:
    """
    parser.add_argument('-qname', help='query species name in orthofinder (see column names of  <Orthogroups.tsv>)')
    parser.add_argument('-qt', help='query species taxid (e.g. use <orthomap qlin -h> to get taxid)')
    parser.add_argument('-sl', help='species list as <orthofinder name><tab><species taxid> '
                                    '(only samples in this list will be processed)')
    parser.add_argument('-oc', help='specify orthofinder <Orthogroups.GeneCounts.tsv> (see Orthogroups directory)')
    parser.add_argument('-og', help='specify orthofinder <Orthogroups.tsv> (see Orthogroups directory)')
    parser.add_argument('-out', help='specify output file <orthomap.tsv> (default: orthomap.tsv)',
                        default='orthomap.tsv')


def get_youngest_common(ql, tl):
    """

    :param ql:
    :param tl:
    :return:
    """
    return [x for x in tl if x in ql][-1]


def get_oldest_common(ql, tl):
    """

    :param ql:
    :param tl:
    :return:
    """
    return ql[min([x for x, y in enumerate(ql) if y in tl])]


def get_orthomap(qname, qt, sl, oc, og, out=None):
    """

    :param qname:
    :param qt:
    :param sl:
    :param oc:
    :param og:
    :param out:
    :return:
    """
    ncbi = NCBITaxa()
    query_lineage = ncbi.get_lineage(qt)
    query_lineage_names_dict = ncbi.get_taxid_translator(query_lineage)
    query_lineage_names = pd.DataFrame([(x, y, query_lineage_names_dict[y]) for x, y in enumerate(query_lineage)])
    query_lineage_names.columns = ['PSnum', 'PStaxID', 'PSname']
    species_list = pd.read_csv(sl, sep='\t', header=None)
    species_list.columns = ['species', 'taxID']
    species_list['lineage'] = species_list.apply(lambda x: ncbi.get_lineage(x[1]), axis=1)
    species_list['youngest_common'] = [get_youngest_common(query_lineage, x) for x in species_list.lineage]
    species_list['youngest_name'] = [list(x.values())[0] for x in [ncbi.get_taxid_translator([x])
                                                                   for x in list(species_list.youngest_common)]]
    print(qname)
    print(qt)
    print(species_list)
    oc_og_dict = {}
    with open(oc, 'r') as oc_lines:
        oc_species = next(oc_lines).strip().split('\t')
        oc_qidx = [x for x, y in enumerate(oc_species) if y == qname]
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
                    oc_og_oldest_common = get_oldest_common(query_lineage, oc_og_hits_youngest_common)
                    oc_og_dict[oc_og[0]] = oc_og_oldest_common
    omap = []
    if out:
        outhandle = open(out, 'w')
        outhandle.write('gene\tOrthogroup\tPSnum\tPStaxID\tPSname\n')
    with open(og, 'r') as og_lines:
        og_species = next(og_lines).strip().split('\t')
        og_qidx = [x for x, y in enumerate(og_species) if y == qname]
        if len(oc_qidx) == 0:
            print('\nError <-qname>: query species name not in orthofinder results, please check spelling\n'
                  'e.g. <head -1 Orthogroups.tsv>')
            sys.exit()
        for og_line in og_lines:
            og_og = og_line.strip().split('\t')
            if og_og[0] not in oc_og_dict:
                continue
            else:
                og_ps = query_lineage_names[query_lineage_names['PStaxID'] ==
                                            oc_og_dict[og_og[0]]].values.tolist()[0]
                og_ps_join = '\t'.join([str(x) for x in og_ps])
                if out:
                    [outhandle.write(x.replace(' ', '') + '\t' + og_og[0] + '\t' + og_ps_join + '\n')
                     for x in og_og[og_qidx[0]].split(',')]
            omap += [[x.replace(' ', ''), og_og[0], og_ps[0], og_ps[1], og_ps[2]]
                     for x in og_og[og_qidx[0]].split(',')]
    if out:
        outhandle.close()
    omap_df = pd.DataFrame(omap)
    omap_df.columns = ['gene', 'Orthogroup', 'PSnum', 'PStaxID', 'PSname']
    return omap_df


def main():
    """

    :return:
    """
    parser = define_parser()
    args = parser.parse_args()
    print(args)
    if not args.qname:
        parser.print_help()
        print('\nError <-qname>: Please specify query species name in orthofinder and taxid')
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
