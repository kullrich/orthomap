#!/usr/bin/python
# -*- coding: UTF-8 -*-


"""
Author: Kristian K Ullrich
date: March 2023
email: ullrich@evolbio.mpg.de
License: GPL-3
"""


import sys
import argparse
import pandas as pd
from ete3 import NCBITaxa, Tree


def define_parser():
    """
    A helper function for using `qlin.py` via the terminal.

    :return: An argparse.ArgumentParser.

    :rtype: argparse.ArgumentParser
    """
    qlin_example = '''qlin example:

    # get query lineage to be used with orthomap later on using query species taxID
    # Mus musculus; 10090
    $ qlin -qt 10090

    # using query species name
    $ qlin -q "Mus musculus"
    '''
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

    :param parser: An argparse.ArgumentParser.

    :type parser: argparse.ArgumentParser
    """
    parser.add_argument('-q', help='query species name')
    parser.add_argument('-qt', help='query species taxID')


def get_qlin(ncbi=None, q=None, qt=None, quiet=False):
    """
    This function searches the NCBI taxonomic database for results matching the
    query name or query taxID.

    Note that if the user specifies both the name and the taxID of a species,
    the returning result is based on the taxID.

    :param ncbi: The NCBI taxonomic database.
    :param q: The name of the queried species.
    :param qt: The taxID of the queried species.
    :param quiet: Specify if output should be quiet.
    :return: A list of information for the queried species such as:
    query name, query taxID, query lineage, query lineage dictonary, query lineage zip,
    query lineage names, reverse query lineage, query kingdom

    :type ncbi: ete3.NCBITaxa
    :type q: str
    :type qt: str
    :type quiet: bool
    :rtype: list

    Example
    --------
    >>> from orthomap import qlin
    >>> from ete3 import NCBITaxa
    >>> qlin.get_qlin(q='Danio rerio')
    """
    qtid = None
    qname = None
    qk = None
    if ncbi is None:
        ncbi = NCBITaxa()
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
    if not quiet:
        print('query name: %s' % qname)
        print('query taxID: %s' % str(qtid))
        print('query kingdom: %s' % qk)
        print(
            'query lineage names: \n%s' % str([qlineagenames_dict[x] + '(' + str(x) + ')' for x in qlineage])
        )
        print('query lineage: \n%s' % str(qlineage))
    return [qname, qtid, qlineage, qlineagenames_dict, qlineagezip, qlineagenames, qlineagerev, qk]


def get_lineage_topo(qt):
    """
    This function returns a species lineage as a tree object for a query species given as taxID.

    :param qt: The taxID of the queried species.
    :return: The lineage of the queried species as an ete3.Tree.

    :type qt: str
    :rtype: ete3.Tree

    Example
    --------
    >>> from orthomap import qlin
    >>> lineage_tree = qlin.get_lineage_topo(qt='10090')
    """
    _, _, _, _, _, qlineagenames, _, _ = get_qlin(qt=qt, quiet=True)
    qln = list(qlineagenames[['PSnum', 'PStaxID', 'PSname']].apply(lambda x: '/'.join(x), axis=1))
    qln = [x.replace('(', '_').replace(')', '_').replace(':', '_') for x in qln]
    tree = Tree('(' * len(qln) + ''.join([str(x) + '),' for x in qln[1::][::-1]])+str(qln[0])+');')
    return tree


def get_youngest_common(ql, tl):
    """
    This function returns the lowest common ancestor (LCA) by comparing the lineage information
    of a query and a target species.

    :param ql: Query species lineage information.
    :param tl: Target species lineage information.
    :return: lowest common ancestor (LCA).

    :type ql: list
    :type tl: list
    :rtype: str

    Example
    --------
    >>> from orthomap import qlin
    >>> # get query species taxonomic lineage information
    >>> _, _, query_lineage, _, _, _, _, _ = qlin.get_qlin(q='Caenorhabditis elegans')
    >>> # get target species taxonomic lineage information
    >>> _, _, target_lineage, _, _, _, _, _ = qlin.get_qlin(q='Mus musculus')
    >>> # get youngest common node
    >>> get_youngest_common(query_lineage, target_lineage)
    """
    return [x for x in tl if x in ql][-1]


def get_oldest_common(ql, tl):
    """
    This function returns the oldest common ancestor (OCA) by comparing the lineage information
    of a query species and a target species.

    The target species can also be a list of LCA values to find the oldest among the given LCA.

    :param ql: Query species lineage information.
    :param tl: Target species lineage information.
    :return: oldest common ancestor (OCA).

    :type ql: list
    :type tl: list
    :rtype: str

    Example
    --------
    >>> from orthomap import qlin
    >>> # get query species taxonomic lineage information
    >>> _, _, query_lineage, _, _, _, _, _ = qlin.get_qlin(q='Caenorhabditis elegans')
    >>> # get target species taxonomic lineage information
    >>> _, _, target_lineage, _, _, _, _, _ = qlin.get_qlin(q='Mus musculus')
    >>> # get oldest common node
    >>> get_oldest_common(query_lineage, target_lineage)
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
        print('\nError <-q> <-qt>: Please specify query species name or taxID')
        sys.exit()
    if args.q and args.qt:
        parser.print_help()
        print('\nWarning: Since both query species name and taxID are given taxID is used')
        sys.exit()
    get_qlin(q=args.q, qt=args.qt, quiet=False)


if __name__ == '__main__':
    main()
