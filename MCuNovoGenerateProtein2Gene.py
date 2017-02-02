# -*- coding: utf-8 -*-
"""
Created on Tue Jan 31 07:45:53 2017

@author: k
extract and save protein to gene info from files
"""
import argparse
parser = argparse.ArgumentParser(description='extract and save protein-to-gene information from files')
parser.add_argument('--cufflinks',nargs = 2, default = None, help = 'two filenames required.\n First, fasta file from typical Cufflinks result. The name line looks like: TCONS_00000001 gene=XLOC_000001\n or TCONS_00000001 XLOC_000001. Second, fasta file translated by transdecoder,\nthe name line looks like: TCONS_00000001|m.1 ...')
parser.add_argument('--denovo',nargs = '+', default = None, help = 'one or more translated proteins by transdecoder,\nthe name line looks like comp0_seq10|m.72 ...\n or looks like: TRINITY_DN100004_c0_g1_i1|m.103699 ...')
parser.add_argument('--maker', nargs = '+', default = None, help = 'one or more protein file similar to maker. \n name line looks like: AGAP000002-RA|m.1 ... \n or looks like: AGAP000002-PA ...')
parser.add_argument('--out',nargs = 1, default = None,help = 'out file name. default is protein2gene.txt in current working directory')

#f = parser.parse_args('--cufflinks cuffT cuffP --maker mak1 mak2'.split())
f = parser.parse_args()

def getProtein2geneCuff(f_cuffTranscript, f_cuffProteins):
    '''
    with two files, return a list of tuples, each with protein ID and gene ID
    two filenames required.\n First, fasta file from typical Cufflinks result. The name line looks like: TCONS_00000001 gene=XLOC_000001\n or TCONS_00000001 XLOC_000001. Second, fasta file translated by transdecoder,\nthe name line looks like: TCONS_00000001|m.1 ...
    (TCONS_00000001|m.1, XLOC_000001)
    '''
    fo = open(f_cuffTranscript,'r')
    dc_transcript2gene = {}
    for line in fo:
        if line[0] == '>':
            transcript, gene = line.split()
            transcript = transcript[1:]
            if 'gene=' in gene:
                gene = gene[5:]
            dc_transcript2gene[transcript] = gene
    fo.close()
    fo = open(f_cuffProteins,'r')
    dc_protein2gene = {}
    for line in fo:
        if line[0] == '>':
            protein = line[1:].split()[0]
            transcript = protein.split('|')[0]
            gene = dc_transcript2gene[transcript]
            dc_protein2gene[protein] = gene
    fo.close()
    return list(dc_protein2gene.items())

def getProtein2geneDenovo(f_denovo):
    '''
    return a list of tuples, each with protein ID and gene ID.
    input is fasta with head line like: comp0_seq10|m.72 or TRINITY_DN100004_c0_g1_i1|m.103699
    '''
    fo = open(f_denovo,'r')
    protein2gene = []
    for line in fo:
        if line[0] == '>':
            protein = line[1:].split()[0]
            transcript = protein.split('|')[0]
            gene = transcript.rsplit('_',1)[0]
            protein2gene.append((protein,gene))
    fo.close()
    return protein2gene

def getProtein2geneMaker(f_maker):
    '''
    return a list of tuples, each with protein ID and gene ID.
    input is fasta with head line like: AGAP000002-RA|m.1 ...  or looks like: AGAP000002-PA ...
    '''
    fo = open(f_maker,'r')
    protein2gene = []
    for line in fo:
        if line[0] == '>':
            protein = line[1:].split()[0]
            gene = protein.split('-')[0]
            protein2gene.append((protein,gene))
    fo.close()
    return protein2gene

p2g = []
if f.cufflinks is not None:
    f_cuffTranscript, f_cuffProteins = f.cufflinks
    p2g = p2g + getProtein2geneCuff(f_cuffTranscript, f_cuffProteins)
if f.denovo is not None:
    for f_denovo in f.denovo:
        p2g = p2g + getProtein2geneDenovo(f_denovo)
if f.maker is not None:
    for f_maker in f.maker:
        p2g = p2g + getProtein2geneMaker(f_maker)

if f.out is None:
    outfile = 'protein2gene.txt'
else:
    outfile = f.out[0]
fout = open(outfile,'w')
for ele in p2g:
    fout.write(ele[0]+'\t'+ele[1]+'\n')
fout.close()

