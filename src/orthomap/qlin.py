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
from ete3 import NCBITaxa


def define_parser():
    """

    :return:
    """
    qlin_example = '''qlin example:

    # get query lineage to be used with orthomap later on using query species taxid
    # Mus musculus; 10090
    qlin -qt 10090
    # using query species name
    qlin -q "Mus musculus"
    '''
    parser = argparse.ArgumentParser(prog='qlin', usage='%(prog)s [options] [<arguments>...]',
                                     description='get query lineage based on ncbi taxonomy', epilog=qlin_example,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    add_argparse_args(parser=parser)
    return parser


def add_argparse_args(parser: argparse.ArgumentParser):
    """

    :param parser:
    :return:
    """
    parser.add_argument('-q', help='query species name')
    parser.add_argument('-qt', help='query species taxid')


def get_qtid(ncbi, q=None, qt=None):
    """

    :param ncbi:
    :param q:
    :param qt:
    :return:
    """
    if q and qt:
        taxid2name = ncbi.get_taxid_translator([int(qt)])
        name2taxid = ncbi.get_name_translator([taxid2name[int(qt)]])
        qname = list(name2taxid.keys())[0]
        qtid = list(taxid2name.keys())[0]
    if not q and qt:
        taxid2name = ncbi.get_taxid_translator([int(qt)])
        name2taxid = ncbi.get_name_translator([taxid2name[int(qt)]])
        qname = list(name2taxid.keys())[0]
        qtid = list(taxid2name.keys())[0]
    if q and not qt:
        name2taxid = ncbi.get_name_translator([q])
        taxid2name = ncbi.get_taxid_translator([name2taxid[q][0]])
        qname = list(name2taxid.keys())[0]
        qtid = list(taxid2name.keys())[0]
    qlineage = ncbi.get_lineage(qtid)
    qlineagenames = ncbi.get_taxid_translator(qlineage)
    qlineagezip = [(a, b) for a, b, in zip(qlineage, [qlineagenames[x] for x in qlineage])]
    qlineagerev = qlineage[::-1]
    if qlineage[2] == 2:
        qk = 'Bacteria'
    if qlineage[2] == 2157:
        qk = 'Archea'
    if qlineage[2] == 2759:
        qk = 'Eukaryota'
    return [qname, qtid, qlineage, qlineagenames, qlineagezip, qlineagerev, qk]


def get_qlin(q=None, qt=None):
    """

    :param q:
    :param qt:
    :return:
    """
    ncbi = NCBITaxa()
    qname, qtid, qlineage, qlineagenames, qlineagezip, qlineagerev, qk = get_qtid(ncbi, q, qt)
    print("query name: %s" % qname)
    print("query taxid: %s" % str(qtid))
    print("query kingdom: %s" % qk)
    print("query lineage names: \n%s" % str([qlineagenames[x] + "(" + str(x) + ")" for x in qlineage]))
    print("query lineage: \n%s" % str(qlineage))
    return [qname, qtid, qlineage, qlineagenames, qlineagezip, qlineagerev, qk]


def main():
    """

    :return:
    """
    parser = define_parser()
    args = parser.parse_args()
    print(args)
    if not args.q and not args.qt:
        parser.print_help()
        print('\nError <-q> <-qt>: Please specify query species name or taxid')
        sys.exit()
    if args.q and args.qt:
        parser.print_help()
        print('\nWarning: Since both query species name and taxid are given taxid is used')
        sys.exit()
    get_qlin(args.q, args.qt)


if __name__ == '__main__':
    main()
