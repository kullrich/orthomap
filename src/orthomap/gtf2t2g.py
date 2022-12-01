#!/usr/bin/python
# -*- coding: UTF-8 -*-


"""
Author: Kristian K Ullrich
date: November 2022
email: ullrich@evolbio.mpg.de
License: GPL-3
"""


import sys
import argparse
import gzip
import pandas as pd


def define_parser():
    """

    :return:
    """
    gtf2t2g_example = '''gtf2t2g example:
    
    # get GTF from Mus musculus
    wget https://ftp.ensembl.org/pub/release-108/gtf/mus_musculus/Mus_musculus.GRCm39.108.chr.gtf.gz
    # create t2g from GTF
    gtf2t2g -i Mus_musculus.GRCm39.108.chr.gtf.gz -o Mus_musculus.GRCm39.108.chr.gtf.t2g.tsv -g -b -p -v -s
    '''
    parser = argparse.ArgumentParser(prog='gtf2t2g', usage='%(prog)s [options] [<arguments>...]',
                                     description='extracts transcript to gene table from GTF', epilog=gtf2t2g_example,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    add_argparse_args(parser=parser)
    return parser


def add_argparse_args(parser: argparse.ArgumentParser):
    """

    :param parser:
    :return:
    """
    parser.add_argument('-i', help='specify GTF input file')
    parser.add_argument('-o', help='specify output file [optional]')
    parser.add_argument('-g', help='specify if gene names should be appended if they exist', action='store_true')
    parser.add_argument('-b', help='specify if gene biotype should be appended if they exist', action='store_true')
    parser.add_argument('-p', help='specify if protein id should be appended if they exist', action='store_true')
    parser.add_argument('-v', help='specify if gene/transcript/protein version should be appended', action='store_true')
    parser.add_argument('-s', help='specify if summary should be printed', action='store_true')
    parser.add_argument('-q', help='specify if output should be quite', action='store_true')


def parse_gtf(gtf, g=False, b=False, p=False, v=False, s=False, output=None, q=False):
    """

    :param gtf: str
        File name of GTF file.
    :param g: bool (default: False)
        Specify if gene names should be appended if they exist.
    :param b: bool (default: False)
        Specify if gene biotype should be appended if they exist.
    :param p: bool (default: False)
        Specify if protein id should be appended if they exist.
    :param v: bool (default: False)
        Specify if gene/transcript/protein version should be appended.
    :param s: bool (default: False)
        Specify if summary should be printed.
    :param output: Optional[str] (default: None)
        File name of output file.
    :param q: bool (default: False)
        Specify if output should be quite.
    :return:

    Example
    --------
    >>> from orthomap import gtf2t2g
    >>> # get gene to transcript table for Danio rerio
    >>> # https://ftp.ensembl.org/pub/release-105/gtf/danio_rerio/Danio_rerio.GRCz11.105.gtf.gz
    >>> query_species_t2g = gtf2t2g.parse_gtf(\
    >>> gtf='Danio_rerio.GRCz11.105.gtf.gz',\
    >>> g=True, b=True, p=True, v=True, s=True, q=True)
    >>> query_species_t2g
    """
    if gtf.endswith('gz'):
        gtf_handle = gzip.open(gtf, 'rt')
    else:
        gtf_handle = open(gtf, 'rt')
    t2g = {}
    t2p = {}
    tc = 0
    gc = 0
    pc = 0
    dc = 0
    for lines in gtf_handle:
        if len(lines) == 0 or lines[0] == '#':
            continue
        line = lines.strip().split('\t')
        if line[2] == 'transcript':
            gid_first = None
            gid_first_version = None
            tid_first = None
            tid_first_version = None
            gname_first = None
            gtype_first = None
            gv_first = None
            tv_first = None
            infosplit = line[8].strip().split(';')
            gid = [x for x in infosplit if 'gene_id' in x]
            if len(gid) > 0:
                if len(gid) == 1:
                    gid_first = gid[0]
                    gid_first = gid_first.replace('gene_id', '').replace(' ', '').replace('"', '')
                else:
                    if not q:
                        print('duplicated gene_id field:\t'+lines)
            else:
                if not q:
                    print('no gene_id field:\t'+lines)
                # continue since no transcript to gene can be evaluated
                continue
            tid = [x for x in infosplit if 'transcript_id' in x]
            if len(tid) > 0:
                if len(tid) == 1:
                    tid_first = tid[0]
                    tid_first = tid_first.replace('transcript_id', '').replace(' ', '').replace('"', '')
                else:
                    if not q:
                        print('duplicated transcript_id field:\t'+lines)
            else:
                if not q:
                    print('no transcript_id field:\t'+lines)
                # continue since no transcript to gene can be evaluated
                continue
            if g:
                gname = [x for x in infosplit if 'gene_name' in x]
                if len(gname) > 0:
                    if len(gname) == 1:
                        gname_first = gname[0]
                        gname_first = gname_first.replace('gene_name', '').replace(' ', '').replace('"', '')
                    else:
                        if not q:
                            print('duplicated gene_name field:\t'+lines)
                else:
                    if not q:
                        print('no gene_name field:\t'+lines)
                    # keep gname_first = None
            if b:
                gtype = [x for x in infosplit if 'gene_biotype' in x]
                if len(gtype) > 0:
                    if len(gtype) == 1:
                        gtype_first = gtype[0]
                        gtype_first = gtype_first.replace('gene_biotype', '').replace(' ', '').replace('"', '')
                    else:
                        if not q:
                            print('duplicated gene_biotype field:\t'+lines)
                else:
                    if not q:
                        print('no gene_biotype field:\t'+lines)
                    # keep gname_gtype = None
            if v:
                gv = [x for x in infosplit if 'gene_version' in x and 'havana_gene_version' not in x]
                if len(gv) > 0:
                    if len(gv) == 1:
                        gv_first = gv[0]
                        gv_first = gv_first.replace('gene_version', '').replace(' ', '').replace('"', '')
                    else:
                        if not q:
                            print('duplicated gene_version field:\t'+lines)
                else:
                    if not q:
                        print('no gene_version field:\t'+lines)
                    # keep gv_first = None
                tv = [x for x in infosplit if 'transcript_version' in x and 'havana_transcript_version' not in x]
                if len(tv) > 0:
                    if len(tv) == 1:
                        tv_first = tv[0]
                        tv_first = tv_first.replace('transcript_version', '').replace(' ', '').replace('"', '')
                    else:
                        if not q:
                            print('duplicated transcript_version field:\t'+lines)
                else:
                    if not q:
                        print('no transcript_version field:\t'+lines)
                    # keep tv_first = None
            if gv_first:
                gid_first_version = gid_first + '.' + gv_first
            if tv_first:
                tid_first_version = tid_first + '.' + tv_first
            if gid_first in t2g:
                if tid_first in t2g[gid_first]:
                    dc += 1
                    if not q:
                        print('duplicated gid-tid: '+gid_first+' '+tid_first)
                    continue
                if tid_first not in t2g[gid_first]:
                    tc += 1
                    t2g[gid_first][tid_first] = [gid_first, str(gid_first_version), tid_first,
                                                 str(tid_first_version), str(gname_first), str(gtype_first)]
            if gid_first not in t2g:
                gc += 1
                tc += 1
                t2g[gid_first] = {}
                t2g[gid_first][tid_first] = [gid_first, str(gid_first_version), tid_first,
                                             str(tid_first_version), str(gname_first), str(gtype_first)]
        if line[2] == 'CDS':
            if p:
                pid_first = None
                pid_first_version = None
                tid_first = None
                pv_first = None
                infosplit = line[8].strip().split(';')
                tid = [x for x in infosplit if 'transcript_id' in x]
                if len(tid) > 0:
                    if len(tid) == 1:
                        tid_first = tid[0]
                        tid_first = tid_first.replace('transcript_id', '').replace(' ', '').replace('"', '')
                    else:
                        if not q:
                            print('duplicated transcript_id field:\t'+lines)
                else:
                    if not q:
                        print('no transcript_id field:\t'+lines)
                    # continue since no transcript to protein can be evaluated
                    continue
                pid = [x for x in infosplit if 'protein_id' in x]
                if len(pid) > 0:
                    if len(pid) == 1:
                        pid_first = pid[0]
                        pid_first = pid_first.replace('protein_id', '').replace(' ', '').replace('"', '')
                    else:
                        if not q:
                            print('duplicated protein_id field:\t'+lines)
                else:
                    if not q:
                        print('no protein_id field:\t'+lines)
                    # continue since no transcript to protein can be evaluated
                    continue
                if v:
                    pv = [x for x in infosplit if 'protein_version' in x]
                    if len(pv) > 0:
                        if len(pv) == 1:
                            pv_first = pv[0]
                            pv_first = pv_first.replace('protein_version', '').replace(' ', '').replace('"', '')
                        else:
                            if not q:
                                print('duplicated protein_version field:\t'+lines)
                    else:
                        if not q:
                            print('no protein_version field:\t'+lines)
                        # keep pv_first = np.nan
                if pv_first:
                    pid_first_version = pid_first + '.' + pv_first
                if tid_first in t2p:
                    continue
                if tid_first not in t2p:
                    pc += 1
                    t2p[tid_first] = [pid_first, str(pid_first_version)]
    for gidk in sorted(t2g.keys()):
        for tidk in sorted(t2g[gidk].keys()):
            pidk = None
            pidk_v = None
            if tidk in t2p:
                pidk = t2p[tidk][0]
                pidk_v = t2p[tidk][1]
            t2g[gidk][tidk] = t2g[gidk][tidk] + [pidk, pidk_v]
            if output:
                output.write('\t'.join([str(x) for x in t2g[gidk][tidk]])+'\n')
    if s:
        print(str(gc)+' gene_id found')
        print(str(tc)+' transcript_id found')
        print(str(tc)+' protein_id found')
        print(str(dc)+' duplicated')
    t2g_df = pd.DataFrame.from_dict({(i, j): t2g[i][j] for i in t2g.keys() for j in t2g[i].keys()},
                                    orient='index', columns=['gene_id', 'gene_id_version',
                                                             'transcript_id', 'transcript_id_version',
                                                             'gene_name', 'gene_type',
                                                             'protein_id', 'protein_id_version'])
    t2g_df.reset_index(drop=True, inplace=True)
    gtf_handle.close()
    return t2g_df


def main():
    """

    :return:
    """
    parser = define_parser()
    args = parser.parse_args()
    print(args)
    if not args.i:
        parser.print_help()
        print('\nError <-i>: Please specify GTF input file')
        sys.exit()
    if args.o:
        output = open(args.o, 'w')
    else:
        output = sys.stdout
    parse_gtf(gtf=args.i, g=args.g, b=args.b, p=args.p, v=args.v,s= args.s, output=output, q=args.q)
    output.close()


if __name__ == '__main__':
    main()
