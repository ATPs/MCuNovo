#!/usr/bin/env python3

import os


description = '''


'''

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-q','--query', help = 'a list of SeqIO or a filename of fasta sequence', required=True)
    parser.add_argument('-s','--subject',help = 'a list of SeqIO or a filename of fasta sequence', required=True)
    parser.add_argument('-m','--min_identity',help='the minimum length of identical AAs required to compare two sequences, default 20', type=int, default=20)
    parser.add_argument('-o', '--outfile', help='whether to put output file. default None, only return a dataframe', default=None)
    parser.add_argument('-e', '--error_rate', help='allowed error rate when counts identical AAs, default 0.02', default=0.02, type=float)
    parser.add_argument('-p', '--muscle', help='location of muscle, default = muscle', default='muscle')
    parser.add_argument('-t', '--threads', help='number of CPUs to use, default=8', default=8, type=int)

    f = parser.parse_args()
    compare2lsSeqs(seqs1=f.query, seqs2=f.subject, identitymin = f.min_identity, outfile = f.outfile, error_rate = f.error_rate, muscle_exe = f.muscle, threads=f.threads)
import argparse






