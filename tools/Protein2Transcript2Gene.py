# -*- coding: utf-8 -*-
import pandas as pd
import re

def getProtein2geneCuff(f_cuffTranscript, f_cuffProteins):
    '''
    two filenames required.
     First, fasta file from typical Cufflinks result. The name line looks like: "TCONS_00000001 gene=XLOC_000001"
     or "TCONS_00000001 XLOC_000001". 
    Second, fasta file translated by transdecoder, the name line looks like: "TCONS_00000001|m.1 ..." or "TCONS_00000001.p5 ..."
    return a dataframe with columns: protein transcript gene
    '''
#    f_cuffTranscript = r"D:\species\genomes\20190622Rhyzopertha_dominica\StringTie\20190625RdoStringTie.fa"
#    f_cuffProteins = r"D:\species\genomes\20190622Rhyzopertha_dominica\StringTie\20190625RdoStringTie.fa.transdecoder.pepUN"
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
    dc_protein2transcript = {}
    for line in fo:
        if line[0] == '>':
            protein = line[1:].split()[0]
            if re.search('\|m\.\d+$', protein):
                transcript = protein.split('|m.')[0]
            else:
                transcript = protein.rsplit('.p', 1)[0]
            dc_protein2transcript[protein] = transcript
    fo.close()
    df = pd.DataFrame()
    df['protein'] = dc_protein2transcript.keys()
    df['transcript'] = dc_protein2transcript.values()
    df['gene'] = df['transcript'].apply(lambda x:dc_transcript2gene[x])
    return df

def getProtein2geneDenovo(f_denovo):
    '''
    one translated proteins by transdecoder,
    the name line looks like: 
        "comp0_seq10|m.72 ..."
        "TRINITY_DN100004_c0_g1_i1|m.103699 ..."
        "TRINITY_DN10002_c0_g1_i1.p1 ..."
        "J21416372.p2 ..."
    return a dataframe with columns: protein transcript gene
    '''
    fo = open(f_denovo,'r')
    proteins = []
    transcripts = []
    genes = []
    for line in fo:
        if line[0] == '>':
            protein = line[1:].split()[0]
            if re.search('\|m\.\d+$', protein):
                transcript = protein.split('|m.')[0]
            else:
                transcript = protein.rsplit('.p', 1)[0]
            gene = transcript.rsplit('_',1)[0]
            proteins.append(protein)
            transcripts.append(transcript)
            genes.append(gene)
    fo.close()
    df = pd.DataFrame()
    df['protein'] = proteins
    df['transcript'] = transcripts
    df['gene'] = genes
    return df

def getProtein2geneMaker(f_maker):
    '''
    one protein file similar to maker. 
     name line looks like: 
         AGAP000002-RA|m.1 ... 
         AGAP000002-RA.p1
         AGAP000002-PA ...
         maker-Scaffold_100_HRSCAF_276-snap-gene-97.64-mRNA-1 ...
    '''
    fo = open(f_maker,'r')
    proteins = []
    transcripts = []
    genes = []
    for line in fo:
        if line[0] == '>':
            protein = line[1:].split()[0]
            if re.search('\|m\.\d+$', protein):
                transcript = protein.split('|m.')[0]
            elif re.search('\.p\d+$', protein):
                transcript = protein.rsplit('.p', 1)[0]
            else:
                transcript = protein
            gene = transcript.rsplit('-', 1)[0]
            proteins.append(protein)
            transcripts.append(transcript)
            genes.append(gene)
    fo.close()
    df = pd.DataFrame()
    df['protein'] = proteins
    df['transcript'] = transcripts
    df['gene'] = genes
    return df

description = '''
extract and save protein-to-gene information from files
return a dataframe with three columns: protein, transcript, and gene
and write a file if --out is provided
'''
des_cuff = '''
works with cufflinks and StringTie output
two filenames required.
 First, fasta file from typical Cufflinks result. The name line looks like: "TCONS_00000001 gene=XLOC_000001"
 or "TCONS_00000001 XLOC_000001". 
Second, fasta file translated by transdecoder, the name line looks like: "TCONS_00000001|m.1 ..." or "TCONS_00000001.p5 ..."
'''
des_denovo = '''
one or more translated proteins by transdecoder,
the name line looks like: 
    "comp0_seq10|m.72 ..."
    "TRINITY_DN100004_c0_g1_i1|m.103699 ..."
    "TRINITY_DN10002_c0_g1_i1.p1 ..."
    "J21416372.p2 ..."
'''
des_maker = '''
one or more protein file similar to maker. 
 name line looks like: 
     AGAP000002-RA|m.1 ... 
     AGAP000002-PA ...
'''

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-c','--cufflinks',nargs = 2, default = None, help = des_cuff)
    parser.add_argument('-d','--denovo',nargs = '+', default = None, help = des_denovo)
    parser.add_argument('-m','--maker', nargs = '+', default = None, help = des_maker)
    parser.add_argument('-o','--out', default = None,help = 'out file name. default is protein2gene.txt in current working directory')
    
    #f = parser.parse_args('--cufflinks cuffT cuffP --maker mak1 mak2'.split())
    f = parser.parse_args()
    
    results = []
    if f.cufflinks is not None:
        f_cuffTranscript, f_cuffProteins = f.cufflinks
        results.append(getProtein2geneCuff(f_cuffTranscript, f_cuffProteins))
    if f.denovo is not None:
        for f_denovo in f.denovo:
            results.append(getProtein2geneDenovo(f_denovo))
    if f.maker is not None:
        for f_maker in f.maker:
            results.append(getProtein2geneMaker(f_maker))
    
    if f.out is None:
        outfile = 'protein2gene.txt'
    else:
        outfile = f.out
    
    df_all = pd.concat(results)
    df_all.to_csv(outfile, sep='\t', index=None)

