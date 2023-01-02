#!/usr/bin/python
# -*- coding: UTF-8 -*-


"""
Author: Kristian K Ullrich
date: November 2022
email: ullrich@evolbio.mpg.de
License: GPL-3
"""


import argparse
import sys

import pandas as pd
from ete3 import NCBITaxa, Tree


def define_parser():
    """
    A helper function for using `qlin.py` via the terminal.

    :return: argparse.ArgumentParser
    """
    qlin_example = """qlin example:

    # get query lineage to be used with orthomap later on using query species taxid
    # Mus musculus; 10090
    $ qlin -qt 10090
    # using query species name
    $ qlin -q "Mus musculus"
    """
    parser = argparse.ArgumentParser(
        prog='qlin',
        usage='%(prog)s [options] [<arguments>...]',
        description='get query lineage based on ncbi taxonomy',
        epilog=qlin_example,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    add_argparse_args(parser=parser)
    return parser


def add_argparse_args(parser: argparse.ArgumentParser):
    """
    This function attaches individual argument specifications to the parser.

    :param parser: argparse.ArgumentParser
    """
    parser.add_argument('-q', help='query species name')
    parser.add_argument('-qt', help='query species taxid')


def get_qtid(ncbi, q=None, qt=None):
    """
    This function searches the NCBI database for results matching the
    query.

    Note that if the user specifies both the name and the taxid of a species,
    the returning result is based on the taxid.

    :param ncbi: ete3.NCBITaxa
        The NCBI database.
    :param q: string
        The name of the queried species.
    :param qt: string
        The taxid of the queried species.
    :return: list
        A list of information for the queried species such as the lineage,
        the lineage names and the kingdom.

    Example
    --------
    >>> >>> from orthomap import qlin
    >>> qlin.get_qlin(q='Danio rerio')
    """
    if qt:
        taxid2name = ncbi.get_taxid_translator([int(qt)])
        qtid, qname = list(taxid2name.items())[0]
    if q and not qt:
        name2taxid = ncbi.get_name_translator([q])
        qname, qtid = list(name2taxid.items())[0]
        qtid = qtid[0]
    qlineage = ncbi.get_lineage(qtid)
    qlineagenames_dict = ncbi.get_taxid_translator(qlineage)
    qlineagezip = [(a, qlineagenames_dict[a]) for a in qlineage]
    qlineagenames = pd.DataFrame([(x, y, qlineagenames_dict[y]) for x, y in enumerate(qlineage)],
                                 columns=['PSnum', 'PStaxID', 'PSname'])
    qlineagenames['PSnum'] = [str(x) for x in list(qlineagenames['PSnum'])]
    qlineagenames['PStaxID'] = [str(x) for x in list(qlineagenames['PStaxID'])]
    qlineagenames['PSname'] = [str(x) for x in list(qlineagenames['PSname'])]
    qlineagerev = qlineage[::-1]
    if qlineage[2] == 2:
        qk = 'Bacteria'
    if qlineage[2] == 2157:
        qk = 'Archea'
    if qlineage[2] == 2759:
        qk = 'Eukaryota'
    return [qname, qtid, qlineage, qlineagenames_dict, qlineagezip, qlineagenames, qlineagerev, qk]


def get_qlin(q=None, qt=None, quite=False):
    """
    A function to retrieve a species' lineage information from the NCBI
    taxonomy database.

    :param q: string
        The name of the queried species.
    :param qt: string
        The taxid of the queried species.
    :param quite: boolean
        Print taxonomic summary.
    :return: list
        A list of information for the queried species such as the lineage,
        the lineage names and the kingdom.

    Example
    --------
    >>> from orthomap import qlin
    >>> # get query species taxonomic lineage information
    >>> query_lineage = qlin.get_qlin(q='Caenorhabditis elegans')
    """
    ncbi = NCBITaxa()
    qname, qtid, qlineage, qlineagenames_dict, qlineagezip, qlineagenames, qlineagerev, qk = get_qtid(ncbi, q, qt)
    if not quite:
        print('query name: %s' % qname)
        print('query taxid: %s' % str(qtid))
        print('query kingdom: %s' % qk)
        print(
            'query lineage names: \n%s' % str([qlineagenames_dict[x] + '(' + str(x) + ')' for x in qlineage])
        )
        print('query lineage: \n%s' % str(qlineage))
    return [qname, qtid, qlineage, qlineagenames_dict, qlineagezip, qlineagenames, qlineagerev, qk]


def get_lineage_topo(qt):
    """
    A function that returns a species lineage as a tree object.

    :param qt: string
        The taxid of the queried species.
    :return: ete3.Tree
        The linage of a species as a ete3.Tree.

    Example
    --------
    >>> from orthomap import qlin
    >>> lineage_tree = qlin.get_lineage_topo(qt=)
    """
    _, _, _, _, _, qlineagenames, _, _ = get_qlin(qt=qt, quite=True)
    qln = list(qlineagenames[['PSnum', 'PStaxID', 'PSname']].apply(lambda x: '/'.join(x), axis=1))
    qln = [x.replace('(', '_').replace(')', '_').replace(':', '_') for x in qln]
    tree = Tree('(' * len(qln) + ''.join([str(x) + '),' for x in qln[1::][::-1]])+str(qln[0])+');')
    return tree


def get_youngest_common(ql, tl):
    """

    :param ql:
    :param tl:
    :return:

    Example
    --------
    >>> from orthomap import qlin
    """
    return [x for x in tl if x in ql][-1]


def get_oldest_common(ql, tl):
    """

    :param ql: 
    :param tl:
    :return:

    Example
    --------
    >>> from orthomap import qlin
    """
    return ql[min([x for x, y in enumerate(ql) if y in tl])]


def main():
    """
    The main function that is being called when `qlin` is used via the terminal.
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
    get_qlin(q=args.q, qt=args.qt, quite=False)


if __name__ == '__main__':
    main()
