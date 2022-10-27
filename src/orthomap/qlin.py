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
    A helper function for using `qlin.py` via the terminal.

    :return: argparse.ArgumentParser
    """
    qlin_example = """qlin example:

    # get query lineage to be used with orthomap later on using query species taxid
    # Mus musculus; 10090
    qlin -qt 10090
    # using query species name
    qlin -q "Mus musculus"
    """
    parser = argparse.ArgumentParser(
        prog="qlin",
        usage="%(prog)s [options] [<arguments>...]",
        description="get query lineage based on ncbi taxonomy",
        epilog=qlin_example,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    add_argparse_args(parser=parser)
    return parser


def add_argparse_args(parser: argparse.ArgumentParser):
    """
    This function attaches individual argument specifications to the parser.

    :param parser: argparse.ArgumentParser
    """
    parser.add_argument("-q", help="query species name")
    parser.add_argument("-qt", help="query species taxid")


def get_qtid(ncbi, q=None, qt=None):
    """
    This function searches the NCBI database for results matching the
    query.

    :param ncbi: ete3.NCBITaxa
        The NCBI database
    :param q: string
        The name of the queried species.
    :param qt: string
        The taxid of the queried species.
    :return: list
        A list of information for the queried species such as the lineage,
        the lineage names and the kingdom.
    """
    # lines 61 (onwards) and 66 (onwards) were the same, because if both arguments
    # are given we are still returning results based on the `taxid`.
    if qt:
        taxid2name = ncbi.get_taxid_translator([int(qt)])
        qtid, qname = list(taxid2name.items())[0]
    if q and not qt:
        name2taxid = ncbi.get_name_translator([q])
        qname, qtid = list(name2taxid.items())[0]
        qtid = qtid[0]
    qlineage = ncbi.get_lineage(qtid)
    qlineagenames = ncbi.get_taxid_translator(qlineage)
    qlineagezip = [(a, qlineagenames[a]) for a in qlineage]  # ToDo verify that this is correct.
    qlineagerev = qlineage[::-1]
    if qlineage[2] == 2:
        qk = "Bacteria"
    if qlineage[2] == 2157:
        qk = "Archea"
    if qlineage[2] == 2759:
        qk = "Eukaryota"
    return [qname, qtid, qlineage, qlineagenames, qlineagezip, qlineagerev, qk]


def get_qlin(q=None, qt=None):
    """
    A function to retrieve a species' lineage information from the NCBI
    taxonomy database.

    :param q: string
        The name of the queried species.
    :param qt: string
        The taxid of the queried species.
    :return: list
        A list of information for the queried species such as the lineage,
        the lineage names and the kingdom.
    """
    ncbi = NCBITaxa()
    qname, qtid, qlineage, qlineagenames, qlineagezip, qlineagerev, qk = get_qtid(ncbi, q, qt)
    print("query name: %s" % qname)
    print("query taxid: %s" % str(qtid))
    print("query kingdom: %s" % qk)
    print(
        "query lineage names: \n%s" % str([qlineagenames[x] + "(" + str(x) + ")" for x in qlineage])
    )
    print("query lineage: \n%s" % str(qlineage))
    return [qname, qtid, qlineage, qlineagenames, qlineagezip, qlineagerev, qk]


def main():
    """
    The main function that is being called when `qlin` is used via the terminal.
    """
    parser = define_parser()
    args = parser.parse_args()
    print(args)
    if not args.q and not args.qt:
        parser.print_help()
        print("\nError <-q> <-qt>: Please specify query species name or taxid")
        sys.exit()
    if args.q and args.qt:
        parser.print_help()
        print("\nWarning: Since both query species name and taxid are given taxid is used")
        sys.exit()
    get_qlin(args.q, args.qt)


if __name__ == "__main__":
    main()
