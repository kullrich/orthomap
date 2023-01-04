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

from ete3 import NCBITaxa


def define_parser():
    """
    A helper function for using `ncbitax.py` via the terminal.

    :return: argparse.ArgumentParser
    """
    ncbitax_example = '''ncbitax example:

    #update ncbi taxonomy database
    ncbitax
    '''
    parser = argparse.ArgumentParser(
        prog='ncbitax',
        usage='%(prog)s [options] [<arguments>...]',
        description='update local ncbi taxonomy database', epilog=ncbitax_example,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    add_argparse_args(parser=parser)
    return parser


def add_argparse_args(parser: argparse.ArgumentParser):
    """
    This function attaches individual argument specifications to the parser.

    :param parser: argparse.ArgumentParser
    """
    parser.add_argument('-u', help='update', action='store_true')


def update_ncbi():
    """
    This function updates or downloads the NCBI taxonomy database using
    the package `ete3`. A parsed version of it will be stored at the home
    directory: `~/.etetoolkit/taxa.sqlite`.

    Example
    --------
    >>> from orthomap import ncbitax
    >>> ncbitax.update_ncbi()
    """
    ncbi = NCBITaxa()
    ncbi.update_taxonomy_database()


def main():
    """
    The main function that is being called when `ncbitax.py` is used via the terminal.
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
