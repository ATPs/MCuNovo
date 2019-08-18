# -*- coding: utf-8 -*-
import pandas as pd
from Bio import SeqIO

def proteinXfinder(filenames, outfile=None):
    '''
    filenames is a list of fasta files, return a dataframe, with protein and count of X
    if outfile is not None, write result to outfile
    '''
    proteins = []
    Xcounts = []
    for filename in filenames:
        for s in SeqIO.parse(filename, 'fasta'):
            seq = str(s.seq)
            seq = seq.upper()
            if 'X' in seq:
                proteins.append(s.id)
                Xcounts.append(seq.count('X'))
    df = pd.DataFrame()
    df['protein'] = proteins
    df['Xcount'] = Xcounts
    if outfile is not None:
        df.to_csv(outfile, sep='\t',index=None)
    return df

description = '''
    filenames is a list of fasta files, return a dataframe, with protein and count of X.
    if outfile is not None, write result to outfile
'''

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-i','--input',nargs = '+', default = None, help = 'files to check', required=True)
    parser.add_argument('-o','--outfile', default = None,help = 'out file name. default is None, only return dataframe')
    f = parser.parse_args()
    proteinXfinder(filenames = f.input, outfile = f.outfile)

