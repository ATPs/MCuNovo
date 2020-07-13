# -*- coding: utf-8 -*-

import pandas as pd
from Bio import SeqIO

def mergeIntervals(intervals):
    '''
    intervals is a list of intervals like [[0, 3], [3, 9], [11, 19], [18, 20], [22, 23]]
    merge overlapped intervals, return [[0,9], [11,20], [22, 23]]
    '''
    sorted_intervals = sorted(intervals) #sort intervals based on start and end
    if len(sorted_intervals) < 2:
        return sorted_intervals
    res = [sorted_intervals[0]]
    for i in range(1, len(sorted_intervals)):
        check = sorted_intervals[i]
        if check[0]<= res[-1][1]:
            if res[-1][1] < check[1]:
                res[-1][1] = check[1]
        else:
            res.append(check)
    return res


def countIntervalLen(intervals):
    '''
    merge intervals and count length. Note, here is for blast6, so [1,3] the length is 3
    '''
    merged_intervals = mergeIntervals(intervals)
    l = 0
    for i in merged_intervals:
        l += (i[1] - i[0] + 1)
    return l


def getML(df):
    '''
    df is a dataframe of blast6, return matched length based on q_start, q_end, s_start s_end
    combine overlapping regions for query and subject
    return matched length in query or subject, which ever is smaller
    '''
    query_intervals = [list(e) for e in zip(df['q_start'], df['q_end'])]
    query_len = countIntervalLen(query_intervals)
    
    subject_intervals =  [list(e) for e in zip(df['s_start'], df['s_end'])]
    subject_len = countIntervalLen(subject_intervals)
    
    identity = (df['identity'] * df['alignment_length']).sum() / df['alignment_length'].sum()
    
    return min(query_len, subject_len), identity

def getSeqLen(seq_ids, f_seq, remove_star = True):
    '''
    f_seq is a filename
    seq_ids is a list/set of ids to work with
    return a dictionary, with sequence ids in seq_id as key, sequence length as value
    if remove_star, the * symbols which indicate end of sequences do not count as part of sequence
    '''
    all_ids = set(seq_ids)
    dc_seqlen = {}
    if f_seq.endswith('.gz'):
        import gzip
        f_seq = gzip.open(f_seq,'rt')
    for s in SeqIO.parse(f_seq, 'fasta'):
        if s.id in all_ids:
            seq = str(s.seq)
            if remove_star:
                if seq[-1] == '*':
                    dc_seqlen[s.id] = len(seq) - 1
                else:
                    dc_seqlen[s.id] = len(seq)
            else:
                dc_seqlen[s.id] = len(seq)
    return dc_seqlen


#f_query = r"D:\species\genomes\20190622Rhyzopertha_dominica\20190622Rhyzopertha_dominica\Rdo_Scaffolds.all.maker.proteins.fasta"
#f_subject = r"D:\species\uniprot\arthropoda\uniprotkb_arthropodaclean"
#f_blast6 = r"D:\species\genomes\20190622Rhyzopertha_dominica\20190714MCuNovo\Maker2U"


def parseBlast6(f_blast6, f_query, f_subject, min_identity=90, outfile=None):
    '''
    Given a filename of Blast6, and query_fasta, subject_fasta file, return a dataframe with columns
    seq1_id seq2_id seq1_len seq2_len match_len
    The content in blast6 looks like:
    AGAP000002-PA	TR:Q7QEI4_ANOGA	100.00	1013	0	0	1	1013	1	1013	0.0	 2086
    which includes the query_id, subject_id, identity, alignment_length, mismatches, gap_opens, q_start, q_end, s_start, s_end, e_value, bit_score.
    match_len is the matching AAs in query or subject whichever smaller after filtering min_identity
    '''
    df_blast = pd.read_csv(f_blast6, sep='\t', header=None)
    df_blast.columns = 'query_id subject_id identity alignment_length mismatches gap_opens q_start q_end s_start s_end e_value bit_score'.split()
    df_match = df_blast.groupby(by=['query_id','subject_id']).apply(getML)
    df_match = df_match.reset_index()
    df_match.columns = 'seq1_id seq2_id temp'.split()
    df_match['match_len'] = df_match['temp'].str[0]
    df_match['identity'] = df_match['temp'].str[1]
    seq1_ids = set(df_match['seq1_id'])
    seq2_ids = set(df_match['seq2_id'])
    dc_seq1len = getSeqLen(seq1_ids, f_query)
    dc_seq2len = getSeqLen(seq2_ids, f_subject)
    df_match['seq1_len'] = df_match['seq1_id'].apply(lambda x:dc_seq1len[x])
    df_match['seq2_len'] = df_match['seq2_id'].apply(lambda x:dc_seq2len[x])

    df_match = df_match['seq1_id seq2_id seq1_len seq2_len match_len identity'.split()]
    
    if outfile is not None:
        df_match.to_csv(outfile, sep='\t', index=None)
    
    return df_match

description = '''    qyery, subject is lists of SeqIO sequences, or fasta filenames
    identitymin is the minimum length of identical AAs required to compare two sequences
    outfile is where to store the result, if None, only return the dataframe
    error_rate is the allowed error rate when counts identical AAs
    muscle_exe is the location of muslce program
    threads is the number of threads to use
    use the output of getPairsFromTwoListSeqs, further calculate matched length
    return a dataframe with five columns: seq1_id, seq2_id, seq1_len, seq2_len, matched_length.
'''

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-i', '--input', help='input blast6 filename')
    parser.add_argument('-q','--query', help = 'filename of query fasta sequence', required=True)
    parser.add_argument('-s','--subject',help = 'filename of subject fasta sequence', required=True)
    parser.add_argument('-m','--min_identity',help='the minimum identity to filter blast6 file , default 95', type=int, default=95)
    parser.add_argument('-o', '--outfile', help='whether to put output file. default None, only return a dataframe', default=None)

    f = parser.parse_args()
    parseBlast6(f_blast6=f.input, f_query=f.query, f_subject=f.subject, min_identity=f.min_identity, outfile=f.outfile)
