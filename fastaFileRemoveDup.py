# -*- coding: utf-8 -*-
"""
Created on Wed Feb  1 13:57:42 2017

@author: k
"""

import MCuNovoGeneSelectorMain

import argparse
parser = argparse.ArgumentParser(description='input a fasta file, output a fasta file with no duplicated sequences')
parser.add_argument('filename',help = 'input file name of fasta file')
parser.add_argument('--error',nargs = 1, default = None, help = 'error rate allowed to determine duplicated sequences, default 0.02')
parser.add_argument('--out',nargs = 1, default = None, help = 'outfile name. Default, input filename +UN')

f = parser.parse_args()
if f.error is None:
    error_rate = 0.02
else:
    error_rate = f.error[0]
if f.out is None:
    outname = f.filename+'UN'
else:
    outname = f.out[0]
f_input = f.filename

ls = MCuNovoGeneSelectorMain.openfile2lsFasta(f_input,'fasta')
lsu = MCuNovoGeneSelectorMain.fasta_within_seq_big_withError(ls,error_rate,6)
fout = open(outname,'w')
for ele in lsu:
    fout.write('>'+ele.description+'\n'+str(ele.seq)+'\n')
fout.close()

