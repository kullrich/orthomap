#!/usr/bin/python
# -*- coding: UTF-8 -*-


"""
Author: Kristian K Ullrich
date: October 2022
email: ullrich@evolbio.mpg.de
License: GPL-3
"""


import argparse
from orthomap import orthomap
from orthomap import qlin
from orthomap import ncbitax


def define_parser():
    parser = argparse.ArgumentParser(prog='orthomap', usage='%(prog)s <sub-command>',
                                     description='orthomap')
    subparsers = parser.add_subparsers(title='sub-commands', help='sub-commands help')
    orthomap_parser = subparsers.add_parser(name='orthomap',
                                            help='extract orthomap from orthofinder output for query species'
                                                 '<orthomap -h>')
    ncbitax_parser = subparsers.add_parser(name='ncbitax', help='update local ncbi taxonomy database <ncbitax -h>')
    qlin_parser = subparsers.add_parser(name='qlin', help='get query lineage based on ncbi taxonomy <qlin -h>')
    orthomap_parser.set_defaults(subcommand='orthomap')
    ncbitax_parser.set_defaults(subcommand='ncbitax')
    qlin_parser.set_defaults(subcommand='qlin')
    orthomap.add_argparse_args(parser=orthomap_parser)
    qlin.add_argparse_args(parser=qlin_parser)
    ncbitax.add_argparse_args(parser=ncbitax_parser)
    return parser


def main():
    parser = define_parser()
    args = parser.parse_args()


if __name__ == '__main__':
    main()
