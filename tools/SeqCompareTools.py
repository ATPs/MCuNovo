# -*- coding: utf-8 -*-
"""
Created on Tue Apr  4 15:17:21 2017
try to work with the xml output of blast directly, to avoid using the BLAST program with edited scoring matrix
@author: ATPs
"""



from Bio.Align.Applications import MuscleCommandline
from Bio import SeqIO
import io
from collections import Counter
import time
import re
import gzip




def openfile2lsFasta(filename,fmt = "fasta", removeStar = True):
    """
    given a file name, return a list of fasta in SeqIO format
    """
    if filename.endswith('.gz'):
        l = list(SeqIO.parse(gzip.open(filename,'rt'),fmt))
    else:
        l = list(SeqIO.parse(open(filename),fmt))
    if removeStar:
        for s in l:
            if s.seq[-1] == '*':
                s.seq = s.seq[:-1]
    return l

def muscleAlignment(seqs, muscle_exe = r"C:\P\Muscle\muscle3.8.31_i86win32.exe"):
    '''
    align sequences with muscle
    given a list of seqs in SeqIO format,
    return a aligned seqs in SeqIO format
    '''
    f_mem = io.StringIO()
    SeqIO.write(seqs,f_mem,'fasta')
    data = f_mem.getvalue()
    muscle_cline = MuscleCommandline(muscle_exe)
    stdout, stderr = muscle_cline(stdin=data)
    return list(SeqIO.parse(io.StringIO(stdout),'fasta'))

def proteinAlignLength(seqs_aln,mincommon = 10,error_rate = 0.02):
    """
    return alignment length based on the output of proteinPairwiseAlignGlobal
    seq1_aln = "MPKSSSN-DLP"
    seq2_aln = "MPRASSNADLP"
    proteinAlignLength(seq1_aln,seq2_aln,1,0), return 8
    proteinAlignLength(seq1_aln,seq2_aln,3,0), return 6
    proteinAlignLength(seq1_aln,seq2_aln,3,0.5), return 10
    """
    seq1_aln = str(seqs_aln[0].seq)
    seq2_aln = str(seqs_aln[1].seq)
    seqlen = len(seq1_aln)
    if seqlen != len(seq2_aln):
        print("input two seqs are not the same length!")
        return None
    seqcommon = "" #find common elements. if two aa not the same, use #. "MP##SSN-AP"
    for num in range(seqlen):
        if seq1_aln[num] == seq2_aln[num]:
            seqcommon = seqcommon + seq1_aln[num]
        else:
            if seq1_aln[num] == "-" or seq2_aln[num] == "-":
                seqcommon = seqcommon + "-"
            else:
                seqcommon = seqcommon + "#"
    seqs = re.split('-+|#{2,}',seqcommon) #split at gap region or if two amino acids are not the same
    common = 0
    for seq in seqs:
        seqlen = len(seq)
        maxerror = int(seqlen * error_rate)
        error = seq.count("#") 
        if error > maxerror:
            seqlen = seqlen - error +maxerror
        if seqlen < mincommon:
            seqlen = 0
        common += seqlen
#        print(seq,seqlen)
    return common

def getProteinAlignLength(seqs, mincommon = 10, error_rate = 0.02,muscle_exe = r"C:\P\Muscle\muscle3.8.31_i86win32.exe"):
    '''
    seqs include two sequences in SeqIO format. Return the aligned length by first align sequences 
    with muscleAlignment, then calculate aligned length with proteinAlignLength
    '''
    seqs_aln = muscleAlignment(seqs,muscle_exe = muscle_exe)
    return proteinAlignLength(seqs_aln, mincommon=mincommon,error_rate= error_rate)

def seqs2kmerdic(seqs,kmerlen=20):
    '''
    seqs is a list of seqIO sequences. 
    Return a dictionary, with kmers with of length 20 (kmerlen) as key, 
    value is a list of id of seqs in seqs
    '''
    from collections import defaultdict
    dckmer = defaultdict(set) #kmer dic, kmer with its seqs
    for n, seq in enumerate(seqs):
        seq = str(seq.seq)
        if len(seq) >= kmerlen:
            for i in range(len(seq)+1-kmerlen):
                kmernum = seq[i:i+kmerlen]
                dckmer[kmernum].add(n)
    #change set to list
    for kmernum in dckmer:
        dckmer[kmernum] = list(dckmer[kmernum])
    return dckmer

def getPairsFromTwoListSeqs(seqs1, seqs2, identitymin = 20, max_target = float('inf')):
    '''
    seqs1 and seqs2 are two list of sequences in SeqIO format
    return a list of tuple, seq1_id, seq2_id, here seq1/2_id means the location of seq in the list of seqs1/2
    for each seq1_id, compare with at most max_target sequences in seqs2
    default, compare all sequences in seqs1 and seqs2 pairs which have at least one identical region with length of 20
    '''
    time1 = time.time()
    df = []
    dc2kmer = seqs2kmerdic(seqs2, kmerlen=identitymin)
    time2 = time.time()
    print('%d seconds, making kmer dictionary for seqs2, kmer length is %d'%(time2 - time1, identitymin))
    for n, seq in enumerate(seqs1):
        s = str(seq.seq)
        #get kmers of s
        skmers = set([s[i:i+identitymin] for i in range(len(s) +1 -identitymin)])
        #get ids in seqs2 which shares a identitymin with s
        s_targets = []
        for kmer in skmers:
            if kmer in dc2kmer:
                s_targets += dc2kmer[kmer]
        if len(s_targets) == 0: # skip if seq have no match in seqs2
            continue
        s_targets = Counter(s_targets)
        s_targets = s_targets.most_common()
        for _n, _target in enumerate(s_targets):
            if _n < max_target:
                df.append((n,_target[0]))
            else:
                break
    time4 = time.time()
    print('%d pairs to compare, total time %d'%(len(df),time4-time1))
    return df


