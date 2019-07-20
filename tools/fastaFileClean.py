# -*- coding: utf-8 -*-
"""
Created on Thu Jan 26 14:48:41 2017

@author: k
"""

def cleanFastaFile(filename, keywords,remove_X = True, outfile = None):
    '''
    Given a filename of fasta file, remove sequences with keyword in keywords or with "X" in protein sequences
    if outfile is None, outfile = filename +'clean'
    '''
    from Bio import SeqIO
    if outfile is None:
        outfile = filename + 'clean'
    fout = open(outfile,'w')
    fo = open(filename,'r')
    total_count = 0
    remove_count = 0
    for sequence in SeqIO.parse(fo,'fasta'):
        total_count += 1
        to_remove = False
        for keyword in keywords:
            keyword = keyword.upper()
            description = sequence.description.upper()
            if keyword in description:
                to_remove = True
                break
        if remove_X:
            seq = sequence.seq
            seq = str(seq)
            seq = seq.upper()
            if 'X' in seq:
                to_remove = True
        if to_remove:
            remove_count += 1
        else:
            SeqIO.write(sequence,fout,'fasta')
    fo.close()
    fout.close()
    print('%d out of %d sequences were removed.'%(remove_count,total_count))

##define parameters
###filename is the input reference fasta file.
#filename = "F:\\Insects\\uniprotkb_arthropoda20170117"
###keywords is a list of keyword that we don't want to see in the head line of fasta sequences
#keywords = ['fragment']
###if remove_X is True, remove sequences with 'X". Can be changed to False
#remove_X = True
###In default, outfile = filename +'clean'. change as you need.
#outfile = None
#
##run the function
#cleanFastaFile(filename, keywords,remove_X, outfile)

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='input a fasta file, output a fasta file with low quality sequences removed')
    parser.add_argument('-k','--keywords', default = 'fragment', help = 'keywords we do not want to see in sequence name, separated by ",". default="fragment"')
    parser.add_argument('-o','--out',help = 'outfile name. Default, input filename +.clean', default=None)
    parser.add_argument('-i','--input',help='input fasta file name', required=True)
    parser.add_argument('-x', '--remove_X', help='whether to remove sequences with X, default = True. True or False, T or F', choices=['True','False', 'T', 'F'], default='True')
    f = parser.parse_args()
    if f.remove_X in ['True','T']:
        remove_X = True
    elif f.remove_X in ['False', 'F']:
        remove_X = False
    else:
        print('wrong format for remove_X')
        exit(0)
    cleanFastaFile(filename=f.input, keywords=f.keywords.split(','), remove_X=remove_X)
    
