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
from Bio import SeqIO
from Bio.Data import CodonTable


def define_parser():
    """
    A helper function for using `cds2aa.py` via the terminal.

    :return: An argparse.ArgumentParser.

    :rtype: argparse.ArgumentParser
    """
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
    parser = argparse.ArgumentParser(
        prog='cds2aa',
        usage='%(prog)s [options] [<arguments>...]',
        description='translate CDS to AA and optional retain longest isoform',
        epilog=cds2aa_example,
        formatter_class=argparse.RawDescriptionHelpFormatter, )
    add_argparse_args(parser=parser)
    return parser


def add_argparse_args(parser: argparse.ArgumentParser):
    """
    This function attaches individual argument specifications to the parser.

    :param parser: An argparse.ArgumentParser.

    :type parser: argparse.ArgumentParser
    """
    parser.add_argument('-i', help='specify fasta input file')
    parser.add_argument('-o', help='specify output file [optional]')
    parser.add_argument('-t', help='transtable [default: std]', default='std')
    parser.add_argument('-r', help='specify CDS source to retain longest isoform')


transtable = {'std': CodonTable.CodonTable(forward_table={
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'TCT': 'S',
    'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'TAT': 'Y', 'TAC': 'Y',
    'TGT': 'C', 'TGC': 'C', 'TGG': 'W', 'CTT': 'L', 'CTC': 'L',
    'CTA': 'L', 'CTG': 'L', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P',
    'CCG': 'P', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'ATT': 'I',
    'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'ACT': 'T', 'ACC': 'T',
    'ACA': 'T', 'ACG': 'T', 'AAT': 'N', 'AAC': 'N', 'AAA': 'K',
    'AAG': 'K', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'GCT': 'A',
    'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'GAT': 'D', 'GAC': 'D',
    'GAA': 'E', 'GAG': 'E', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G',
    'GGG': 'G', '---': 'X', '--A': 'X', '--C': 'X', '--G': 'X',
    '--T': 'X', '-A-': 'X', '-AA': 'X', '-AC': 'X', '-AG': 'X',
    '-AT': 'X', '-C-': 'X', '-CA': 'X', '-CC': 'X', '-CG': 'X',
    '-CT': 'X', '-G-': 'X', '-GA': 'X', '-GC': 'X', '-GG': 'X',
    '-GT': 'X', '-T-': 'X', '-TA': 'X', '-TC': 'X', '-TG': 'X',
    '-TT': 'X', 'A--': 'X', 'A-A': 'X', 'A-C': 'X', 'A-G': 'X',
    'A-T': 'X', 'AA-': 'X', 'AC-': 'X', 'AG-': 'X', 'AT-': 'X',
    'C--': 'X', 'C-A': 'X', 'C-C': 'X', 'C-G': 'X', 'C-T': 'X',
    'CA-': 'X', 'CC-': 'X', 'CG-': 'X', 'CT-': 'X', 'G--': 'X',
    'G-A': 'X', 'G-C': 'X', 'G-G': 'X', 'G-T': 'X', 'GA-': 'X',
    'GC-': 'X', 'GG-': 'X', 'GT-': 'X', 'T--': 'X', 'T-A': 'X',
    'T-C': 'X', 'T-G': 'X', 'T-T': 'X', 'TA-': 'X', 'TC-': 'X',
    'TG-': 'X', 'TT-': 'X', 'NNN': 'X', 'GCN': 'A', 'CGN': 'R',
    'MGR': 'R', 'AAY': 'N', 'GAY': 'D', 'TGY': 'C', 'CAR': 'Q',
    'GAR': 'E', 'GGN': 'G', 'CAY': 'H', 'ATH': 'I', 'YTR': 'L',
    'CTN': 'L', 'AAR': 'K', 'TTY': 'F', 'CCN': 'P', 'TCN': 'S',
    'AGY': 'S', 'ACN': 'T', 'TAY': 'Y', 'GTN': 'V', 'TAR': '*',
    'TRA': '*'},
    stop_codons=['TAA', 'TAG', 'TGA', ],
    start_codons=['TTG', 'CTG', 'ATG', ])}


def get_gene(description,
             source):
    """
    This function extracts the gene ID from a sequence description based on a given source.

    :param description: Sequence description.
    :param source: Sequence source.
    :return: GeneID.

    :type description: str
    :type source: str
    :rtype: str
    """
    gene = None
    if source == 'NCBI':
        gene = description.split('[gene=')[1].split(' ')[0].replace(']', '')
    if source == 'ENSEMBL':
        gene = description.split('gene:')[1].split(' ')[0]
    if source == 'WORMBASE':
        gene = description.split('gene=')[1].split(' ')[0]
    return gene


def get_gene_len_dict(record_iter,
                      source):
    """
    This function creates a sequence length dictionary.

    :param record_iter: Sequence iterable.
    :param source: Sequence source.
    :return: Sequence length dictionary.

    :type record_iter: Bio.SeqIO.FastaIO.FastaIterator
    :type source: str
    :rtype: dictionary
    """
    record_dict = {}
    for record in record_iter:
        gene = get_gene(record.description, source)
        if gene in record_dict:
            if record_dict[gene][0] < len(record):
                record_dict[gene] = [len(record), record]
        if gene not in record_dict:
            record_dict[gene] = [len(record), record]
    return record_dict


def cds2aa_record(record_iter,
                  args):
    """
    This function translates nucleotide to amino acids assuming that cds is in frame 0.

    This is a generator function, the records argument should
    be a list or iterator returning SeqRecord objects.

    :param record_iter: Sequence iterable.
    :param args: An argparse.Namespace.
    :return: Sequence object.

    :type record_iter: Bio.SeqIO.FastaIO.FastaIterator
    :type args: argparse.Namespace
    :rtype: Bio.SeqIO.SeqRecord
    """
    for record in record_iter:
        aa = SeqIO.SeqRecord(record.seq.translate(transtable[args.t]), name=record.name, id=record.name,
                             description=record.name)
        yield aa


def cds2aa_fasta(args,
                 parser):
    """
    This function

    :param args:
    :param parser:
    :return:

    :type args:
    :type parser:
    :rtype:
    """
    record_iter = None
    if args.i is None and sys.stdin.isatty():
        parser.print_help()
        sys.exit('\nPlease provide STDIN or input file')
    if args.i is None and not sys.stdin.isatty():
        record_iter = SeqIO.parse(sys.stdin, "fasta")
    else:
        record_iter = SeqIO.parse(args.i, "fasta")
    if args.r:
        record_gene_len_dict = get_gene_len_dict(record_iter, args.r)
        record_iter = iter([x[1] for x in record_gene_len_dict.values()])
    cds2aa_iter = cds2aa_record(record_iter, args)
    if args.o is None:
        SeqIO.write(cds2aa_iter, sys.stdout, "fasta")
    else:
        count = SeqIO.write(cds2aa_iter, args.o, "fasta")
        print("translated %i sequences" % count)


def main():
    """
    The main function that is being called when `cds2aa` is used via the terminal.
    """
    parser = define_parser()
    args = parser.parse_args()
    if args.o is None:
        sys.stderr.write(str(args))
    else:
        print(args)
    cds2aa_fasta(args, parser)


if __name__ == '__main__':
    main()
