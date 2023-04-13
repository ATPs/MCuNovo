# -*- coding: utf-8 -*-
"""
Created on Wed Feb  1 13:57:42 2017

@author: k
"""
from Bio import SeqIO
import time
try:
    import tqdm
    tqdm_exist = True
except:
    tqdm_exist = False

def errorMatch(seq1, seq2, errors=2):
    """
    Given seq1 and seq2, len(seq1) <= len(seq2), return whether they match each other with allowed error.
    """
    if len(seq1) > len(seq2):
        return False
    step = len(seq1)//(errors+1)
    if step == 0:
        return True
    if errors == 0:
        return seq1 in seq2
    parts = [seq1[i:i+step] for i in range(0,len(seq1),step)] #separate seq1 to error+1 parts
    if len(parts[-1]) < step:
        parts[-2] = parts[-2]+parts[-1]
        parts.pop()
    similar = False
    sameslist =[]
    for i in range(errors+1):
        findsame = seq2.find(parts[i])
        if  findsame >= 0:
            similar = True
            sameslist.append((i,findsame))
    if similar == False:
        return False
    for i,j in sameslist:
        if j-step*(i)>=0 and j-step*(i)+len(seq1) <= len(seq2):
            seq2n = seq2[j-step*i:j-step*i +len(seq1)]
            missmatched = 0
            for k in range(len(seq1)):
                if seq1[k] != seq2n[k]:
                    missmatched += 1
                if missmatched > errors:
                    break
            if missmatched <= errors:
                return True
    return False

def getDicKmernum(myfasta, kmerlen = 6):
    '''
    myfasta is a list of SeqIO elements
    return a dictionary, with kmer with length of kmerlen as key, and a list of index of sequences with that kmer
    '''
    dickmernum = {} #kmer dic, kmer with its seqs
    for dummyi in range(len(myfasta)):
        seq = str(myfasta[dummyi].seq)
        for i in range(len(seq)+1-kmerlen):
            kmernum = seq[i:i+kmerlen]
            if kmernum not in dickmernum:
                dickmernum[kmernum] = []
            dickmernum[kmernum].append(dummyi)
    
    time1 = time.time() #change values of dickmernum to list
    for kmernum in dickmernum:
        dickmernum[kmernum] = list(set(dickmernum[kmernum]))
    print(time.time()-time1)
    
    return dickmernum


def fasta_within_seq_big_withError(myfasta, error_rate = 0.02,kmerlen = 6):
    """
    myfasta is a list of SeqIO elements
    if a sequence is part of the other, with error_rate allowed, then remove this sequence.
    return a list of non-redundant SeqIO fasta
    """
    time1 = time.time()
    dickmernum = getDicKmernum(myfasta, kmerlen = kmerlen)
    print(time.time()-time1)
    
    toremove = set()
    from collections import Counter
    if tqdm_exist:
        to_iter = tqdm.tqdm(range(len(myfasta)))
    else:
        to_iter = range(len(myfasta))
    for num1 in to_iter:
        seq1 = str(myfasta[num1].seq)
        seq1kmers = set() # all kmernum, here is kmer5 in seq1
        for i in range(len(seq1)+1-kmerlen):
            seq1kmers.add(seq1[i:i+kmerlen])
    #    print(time.time()-time1)
        seq1targets = []
        for kmernum in seq1kmers:
            seq1targets += dickmernum[kmernum]
        seq1targets = Counter(seq1targets) # count the number of common kmers for each targets
        seq1targets = seq1targets.most_common() # sort the targets based on the number of commn kmers
    #    print(time.time()-time1)
        errors = int(len(seq1)*error_rate)
        for seq2id, seq2_counts in seq1targets:
            if seq2id != num1:
                seq2 = str(myfasta[seq2id].seq)
                if seq2_counts >= len(seq1kmers)-errors*kmerlen:
                    if len(seq1) <= len(seq2):
                        if seq2id not in toremove:
                            if seq2_counts >=2:
                                if errorMatch(seq1,seq2,errors):
                                    toremove.add(num1)
                                    break
    
    print(time.time()-time1)
    print('total removed sequence number is')
    print(len(toremove))
    nonredunfasta =[]
    for i in range(len(myfasta)):
        if i not in toremove:
            nonredunfasta.append(myfasta[i])
    return nonredunfasta

description='''
input a fasta file, output a fasta file with no duplicated sequences
default allowed error rate is 0.02 when align two sequences together.
'''

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('filename',help = 'input file name of fasta file')
    parser.add_argument('-e', '--error', default = None, help = 'error rate allowed to determine duplicated sequences, default 0.02')
    parser.add_argument('-k', '--kmerlen', default = 15, help = 'kmer length used to build a dict for sequences. default: 15', type=int)
    parser.add_argument('-o','--out', default = None, help = 'outfile name. Default, input filename +UN')
    
    f = parser.parse_args()
    if f.error is None:
        error_rate = 0.02
    else:
        error_rate = float(f.error)
    if f.out is None:
        outname = f.filename+'UN'
    else:
        outname = f.out
    f_input = f.filename
    
    ls = list(SeqIO.parse(f_input,'fasta'))
    print('finished reading the fasta file')
    lsu = fasta_within_seq_big_withError(ls,error_rate,f.kmerlen)
    fout = open(outname,'w')
    for ele in lsu:
        fout.write('>'+ele.description+'\n'+str(ele.seq)+'\n')
    fout.close()

