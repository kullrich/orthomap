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
import gzip


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


def parse_gtf(gtf, output, g, b, p, v, s):
    t2g = {}
    t2p = {}
    tc = 0
    gc = 0
    pc = 0
    dc = 0
    for lines in gtf:
        if len(lines) == 0 or lines[0] == '#':
            continue
        line = lines.strip().split('\t')
        if line[2] == 'transcript':
            infosplit = line[8].strip().split(';')
            gid = [x for x in infosplit if 'gene_id' in x]
            if len(gid) > 0:
                if len(gid) == 1:
                    gid_first = gid[0]
                    gid_first = gid_first.replace('gene_id', '').replace(' ', '').replace('"', '')
                else:
                    print('duplicated gene_id field:\t'+lines)
            else:
                print('no gene_id field:\t'+lines)
                continue
            tid = [x for x in infosplit if 'transcript_id' in x]
            if len(tid) > 0:
                if len(tid) == 1:
                    tid_first = tid[0]
                    tid_first = tid_first.replace('transcript_id', '').replace(' ', '').replace('"', '')
                else:
                    print('duplicated transcript_id field:\t'+lines)
            else:
                print('no transcript_id field:\t'+lines)
                continue
            if g:
                gname = [x for x in infosplit if 'gene_name' in x]
                if len(gname) > 0:
                    if len(gname) == 1:
                        gname_first = gname[0]
                        gname_first = gname_first.replace('gene_name', '').replace(' ', '').replace('"', '')
                    else:
                        print('duplicated gene_name field:\t'+lines)
                else:
                    gname_first = ''
            if b:
                gtype = [x for x in infosplit if 'gene_biotype' in x]
                if len(gtype) > 0:
                    if len(gtype) == 1:
                        gtype_first = gtype[0]
                        gtype_first = gtype_first.replace('gene_biotype', '').replace(' ', '').replace('"', '')
                    else:
                        print('duplicated gene_biotype field:\t'+lines)
                else:
                    gtype_first = ''
            if v:
                gv = [x for x in infosplit if 'gene_version' in x and 'havana_gene_version' not in x]
                if len(gv) > 0:
                    if len(gv) == 1:
                        gv_first = gv[0]
                        gv_first = gv_first.replace('gene_version', '').replace(' ', '').replace('"', '')
                    else:
                        print('duplicated gene_version field:\t'+lines)
                else:
                    print('no gene_version field:\t'+lines)
                    continue
                gid_first = gid_first+'.'+gv_first
                tv = [x for x in infosplit if 'transcript_version' in x and 'havana_transcript_version' not in x]
                if len(tv) > 0:
                    if len(tv) == 1:
                        tv_first = tv[0]
                        tv_first = tv_first.replace('transcript_version', '').replace(' ', '').replace('"', '')
                    else:
                        print('duplicated transcript_version field:\t'+lines)
                else:
                    print('no transcript_version field:\t'+lines)
                    continue
                tid_first = tid_first+'.'+tv_first
            if gid_first in t2g:
                if tid_first in t2g[gid_first]:
                    dc += 1
                    print('duplicated gid-tid: '+gid_first+' '+tid_first)
                    continue
                if tid_first not in t2g[gid_first]:
                    tc += 1
                    if g and b:
                        t2g[gid_first][tid_first] = [gid_first, tid_first, gname_first, gtype_first]
                    if g and not b:
                        t2g[gid_first][tid_first] = [gid_first, tid_first, gname_first]
                    if not g and not b:
                        t2g[gid_first][tid_first] = [gid_first, tid_first]
            if gid_first not in t2g:
                gc += 1
                tc += 1
                t2g[gid_first] = {}
                if g and b:
                    t2g[gid_first][tid_first] = [gid_first, tid_first, gname_first, gtype_first]
                if g and not b:
                    t2g[gid_first][tid_first] = [gid_first, tid_first, gname_first]
                if not g and b:
                    t2g[gid_first][tid_first] = [gid_first, tid_first, gtype_first]
                if not g and not b:
                    t2g[gid_first][tid_first] = [gid_first, tid_first]
        if line[2] == 'CDS':
            if p:
                infosplit = line[8].strip().split(';')
                tid = [x for x in infosplit if 'transcript_id' in x]
                if len(tid) > 0:
                    if len(tid) == 1:
                        tid_first = tid[0]
                        tid_first = tid_first.replace('transcript_id', '').replace(' ', '').replace('"', '')
                    else:
                        print('duplicated transcript_id field:\t'+lines)
                else:
                    print('no transcript_id field:\t'+lines)
                    continue
                pid = [x for x in infosplit if 'protein_id' in x]
                if len(pid) > 0:
                    if len(pid) == 1:
                        pid_first = pid[0]
                        pid_first = pid_first.replace('protein_id', '').replace(' ', '').replace('"', '')
                    else:
                        print('duplicated protein_id field:\t'+lines)
                else:
                    print('no protein_id field:\t'+lines)
                    continue
                if v:
                    tv = [x for x in infosplit if 'transcript_version' in x and 'havana_transcript_version' not in x]
                    if len(tv) > 0:
                        if len(tv) == 1:
                            tv_first = tv[0]
                            tv_first = tv_first.replace('transcript_version', '').replace(' ', '').replace('"', '')
                        else:
                            print('duplicated transcript_version field:\t'+lines)
                    else:
                        print('no transcript_version field:\t'+lines)
                        continue
                    tid_first = tid_first+'.'+tv_first
                    pv = [x for x in infosplit if 'protein_version' in x]
                    if len(pv) > 0:
                        if len(pv) == 1:
                            pv_first = pv[0]
                            pv_first = pv_first.replace('protein_version', '').replace(' ', '').replace('"', '')
                        else:
                            print('duplicated protein_version field:\t'+lines)
                    else:
                        print('no protein_version field:\t'+lines)
                        continue
                    pid_first = pid_first+'.'+pv_first
                if tid_first in t2p:
                    continue
                if tid_first not in t2p:
                    pc += 1
                    t2p[tid_first] = pid_first
    if p:
        for gidk in sorted(t2g.keys()):
            for tidk in sorted(t2g[gidk].keys()):
                pidk = ''
                if tidk in t2p:
                    pidk = t2p[tidk]
                output.write('\t'.join(t2g[gidk][tidk])+'\t'+pidk+'\n')
    else:
        for gidk in sorted(t2g.keys()):
            for tidk in sorted(t2g[gidk].keys()):
                output.write('\t'.join(t2g[gidk][tidk])+'\n')
    if s:
        print(str(gc)+' gene_id found')
        print(str(tc)+' transcript_id found')
        print(str(tc)+' protein_id found')
        print(str(dc)+' duplicated')


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
    if args.i.endswith('gz'):
        gtf = gzip.open(args.i, 'rt')
    else:
        gtf = open(args.i, 'rt')
    if args.o:
        output = open(args.o, 'w')
    else:
        output = sys.stdout
    parse_gtf(gtf, output, args.g, args.b, args.p, args.v, args.s)
    gtf.close()
    output.close()


if __name__ == '__main__':
    main()
