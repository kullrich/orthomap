#!/usr/bin/python
# -*- coding: UTF-8 -*-


"""
Author: Kristian K Ullrich
date: January 2023
email: ullrich@evolbio.mpg.de
License: GPL-3
"""


import argparse
from orthomap import of2orthomap, ncbitax, qlin


def define_parser():
    parser = argparse.ArgumentParser(prog='orthomap', usage='%(prog)s <sub-command>',
                                     description='orthomap')
    subparsers = parser.add_subparsers(title='sub-commands', help='sub-commands help')
    of2orthomap_example = '''of2orthomap example:

    #
    of2orthomap -seqname -qt 10090 -sl -oc -og
    '''
    of2orthomap_parser = subparsers.add_parser(name='of2orthomap',
                                               help='extract orthomap from OrthoFinder output for query species'
                                                    '<orthomap -h>',
                                               epilog=of2orthomap_example,
                                               formatter_class=argparse.RawDescriptionHelpFormatter)
    ncbitax_parser = subparsers.add_parser(name='ncbitax', help='update local ncbi taxonomy database <ncbitax -h>')
    qlin_parser = subparsers.add_parser(name='qlin', help='get query lineage based on ncbi taxonomy <qlin -h>')
    of2orthomap_parser.set_defaults(subcommand='of2orthomap')
    ncbitax_parser.set_defaults(subcommand='ncbitax')
    qlin_parser.set_defaults(subcommand='qlin')
    of2orthomap.add_argparse_args(parser=of2orthomap_parser)
    qlin.add_argparse_args(parser=qlin_parser)
    ncbitax.add_argparse_args(parser=ncbitax_parser)
    return parser


def main():
    parser = define_parser()
    args = parser.parse_args()


if __name__ == '__main__':
    main()
