#!/usr/bin/python
# -*- coding: UTF-8 -*-


"""
Author: Kristian K Ullrich
date: March 2023
email: ullrich@evolbio.mpg.de
License: GPL-3
"""


import os
import sys
import argparse
from orthomap import cds2aa, gtf2t2g, ncbitax, of2orthomap, qlin


def define_parser():
    """
    A helper function for using `orthomap` via the terminal.

    :return: An argparse.ArgumentParser.

    :rtype: argparse.ArgumentParser
    """
    parser = argparse.ArgumentParser(prog='orthomap', usage='%(prog)s <sub-command>',
                                     description='orthomap')
    subparsers = parser.add_subparsers(dest='subcommand', title='sub-commands', help='sub-commands help')
    cds2aa_example = '''cds2aa example:

    # to get CDS from Danio rerio on Linux run:
    $ wget https://ftp.ensembl.org/pub/release-105/fasta/danio_rerio/cds/Danio_rerio.GRCz11.cds.all.fa.gz
    $ gunzip Danio_rerio.GRCz11.cds.all.fa.gz

    # on Mac:
    $ curl https://ftp.ensembl.org/pub/release-105/fasta/danio_rerio/cds/Danio_rerio.GRCz11.cds.all.fa.gz --remote-name
    $ gunzip Danio_rerio.GRCz11.cds.all.fa.gz
    
    # translate and retain longest isoform from CDS fasta file:
    $ cds2aa -i Danio_rerio.GRCz11.cds.all.fa -r ENSEMBL -o Danio_rerio.GRCz11.aa.all.longest.fa
    '''
    gtf2t2g_example = '''gtf2t2g example:

    # to get GTF from Mus musculus on Linux run:
    $ wget https://ftp.ensembl.org/pub/release-108/gtf/mus_musculus/Mus_musculus.GRCm39.108.chr.gtf.gz

    # on Mac:
    $ curl https://ftp.ensembl.org/pub/release-108/gtf/mus_musculus/Mus_musculus.GRCm39.108.chr.gtf.gz --remote-name

    # create t2g from GTF:
    $ gtf2t2g -i Mus_musculus.GRCm39.108.chr.gtf.gz -o Mus_musculus.GRCm39.108.chr.gtf.t2g.tsv -g -b -p -v -s
    '''
    ncbitax_example = '''ncbitax example:

    #update ncbi taxonomy database:
    ncbitax -u
    '''
    of2orthomap_example = '''of2orthomap example:

    # download OrthoFinder example:
    $ wget https://github.com/kullrich/orthomap/raw/main/examples/ensembl_105_orthofinder_Orthogroups.GeneCount.tsv.zip
    $ wget https://github.com/kullrich/orthomap/raw/main/examples/ensembl_105_orthofinder_Orthogroups.tsv.zip
    $ wget https://github.com/kullrich/orthomap/raw/main/examples/ensembl_105_orthofinder_species_list.tsv

    # extract orthomap:
    $ of2orthomap -seqname Danio_rerio.GRCz11.cds.longest -qt 7955 \\
      -sl ensembl_105_orthofinder_species_list.tsv \\
      -oc ensembl_105_orthofinder_Orthogroups.GeneCount.tsv.zip \\
      -og ensembl_105_orthofinder_Orthogroups.tsv.zip
    '''
    qlin_example = '''qlin example:

    # get query lineage to be used with orthomap later on using query species taxid
    # Mus musculus; 10090
    $ qlin -qt 10090

    # using query species name
    $ qlin -q "Mus musculus"
    '''
    cds2aa_parser = subparsers.add_parser(name='cds2aa',
                                          help='translate CDS to AA and optional retain longest isoform <cds2aa -h>',
                                          epilog=cds2aa_example,
                                          formatter_class=argparse.RawDescriptionHelpFormatter)
    gtf2t2g_parser = subparsers.add_parser(name='gtf2t2g',
                                           help='extracts transcript to gene table from GTF <gtf2t2g -h>',
                                           epilog=gtf2t2g_example,
                                           formatter_class=argparse.RawDescriptionHelpFormatter)
    ncbitax_parser = subparsers.add_parser(name='ncbitax',
                                           help='update local ncbi taxonomy database <ncbitax -h>',
                                           epilog=ncbitax_example,
                                           formatter_class=argparse.RawDescriptionHelpFormatter)
    of2orthomap_parser = subparsers.add_parser(name='of2orthomap',
                                               help='extract orthomap from OrthoFinder output for query species'
                                                    '<orthomap -h>',
                                               epilog=of2orthomap_example,
                                               formatter_class=argparse.RawDescriptionHelpFormatter)
    qlin_parser = subparsers.add_parser(name='qlin',
                                        help='get query lineage based on ncbi taxonomy <qlin -h>',
                                        epilog=qlin_example,
                                        formatter_class=argparse.RawDescriptionHelpFormatter)
    cds2aa.add_argparse_args(parser=cds2aa_parser)
    gtf2t2g.add_argparse_args(parser=gtf2t2g_parser)
    ncbitax.add_argparse_args(parser=ncbitax_parser)
    of2orthomap.add_argparse_args(parser=of2orthomap_parser)
    qlin.add_argparse_args(parser=qlin_parser)
    return parser


def main():
    """
    The main function that is being called when `orthomap` is used via the terminal.
    """
    parser = define_parser()
    args = parser.parse_args()
    if args.subcommand is None:
        parser.print_help()
        sys.exit()
    print(args)
    if args.subcommand == 'cds2aa':
        if args.o is None:
            sys.stderr.write(str(args))
        else:
            print(args)
        cds2aa.cds2aa_fasta(args, parser)
    if args.subcommand == 'gtf2t2g':
        if not args.i:
            parser.print_help()
            print('\nError <-i>: Please specify GTF input file')
            sys.exit()
        if args.o:
            if os.path.exists(args.o) and not args.overwrite:
                print('\nError <-overwrite>: output file exists, please set to True if it should be overwritten\n')
                sys.exit()
            output = open(args.o, 'w')
        else:
            output = sys.stdout
        gtf2t2g.parse_gtf(gtf=args.i, g=args.g, b=args.b, p=args.p, v=args.v, s=args.s, output=output, q=args.q)
        output.close()
    if args.subcommand == 'ncbitax':
        if not args.u:
            parser.print_help()
            print('\nError <-u>: Please specify if you like to update <-u>')
            sys.exit()
        if args.u:
            ncbitax.update_ncbi()
    if args.subcommand == 'of2orthomap':
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
        of2orthomap.get_orthomap(seqname=args.seqname, qt=args.qt, sl=args.sl, oc=args.oc, og=args.og, out=args.out,
                                 quiet=False, continuity=True, overwrite=args.overwrite)
    if args.subcommand == 'qlin':
        if not args.q and not args.qt:
            parser.print_help()
            print('\nError <-q> <-qt>: Please specify query species name or taxid')
            sys.exit()
        if args.q and args.qt:
            parser.print_help()
            print('\nWarning: Since both query species name and taxid are given taxid is used')
            sys.exit()
        qlin.get_qlin(q=args.q, qt=args.qt, quiet=False)


if __name__ == '__main__':
    main()
