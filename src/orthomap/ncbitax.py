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
    ncbitax_example = '''ncbitax example:

    #update ncbi taxonomy database
    ncbitax
    '''
    parser = argparse.ArgumentParser(prog='ncbitax', usage='%(prog)s [options] [<arguments>...]',
                                     description='update local ncbi taxonomy database', epilog=ncbitax_example,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    add_argparse_args(parser=parser)
    return parser


def add_argparse_args(parser: argparse.ArgumentParser):
    """

    :param parser:
    :return:
    """
    parser.add_argument('-u', help='update', action='store_true')


def update_ncbi():
    """

    :return:
    """
    ncbi = NCBITaxa()
    ncbi.update_taxonomy_database()


def main():
    """
    update local ncbi taxonomy database

    :return:
    """
    parser = define_parser()
    args = parser.parse_args()
    print(args)
    if not args.u:
        parser.print_help()
        print('\nError <-u>: Please specify if you like to update <-u>')
        sys.exit()
    if args.u:
        update_ncbi()


if __name__ == '__main__':
    main()
