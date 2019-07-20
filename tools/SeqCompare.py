from multiprocessing import Pool
import sys
sys.path.append(__file__)
import pandas as pd

from SeqCompareTools import openfile2lsFasta, getProteinAlignLength,getPairsFromTwoListSeqs

def compare2lsSeqs(seqs1, seqs2, identitymin = 20, outfile = None,error_rate = 0.02,muscle_exe = r"C:\P\Muscle\muscle3.8.31_i86win32.exe",threads=12):
    '''
    seqs1, seqs2 is lists of SeqIO sequences, or fasta filenames
    identitymin is the minimum length of identical AAs required to compare two sequences
    outfile is where to store the result, if None, only return the dataframe
    error_rate is the allowed error rate when counts identical AAs
    muscle_exe is the location of muslce program
    threads is the number of threads to use
    use the output of getPairsFromTwoListSeqs, further calculate matched length
    return a dataframe with three columns: seq1_id, seq2_id, seq1_len, seq2_len, matched_length.
    '''
    # convert seqs1, seqs2 to list of SeqIO if they are not
    if isinstance(seqs1,str):
        seqs1 = openfile2lsFasta(seqs1)
    if isinstance(seqs2,str):
        seqs2 = openfile2lsFasta(seqs2)

    
    # sequence pares to compare
    pairs = getPairsFromTwoListSeqs(seqs1, seqs2, identitymin = identitymin, max_target = float('inf'))
    print('total pairs to compare:', len(pairs))
    
    parameters = [[[seqs1[pair[0]], seqs2[pair[1]]], identitymin, error_rate,muscle_exe] for pair in pairs]
    pool = Pool(threads)
    match_lengths = pool.starmap(getProteinAlignLength,parameters)
    pool.close()
    
    df = pd.DataFrame()
    df['seq1_id'] = [seqs1[i[0]].id for i in pairs]
    df['seq2_id'] = [seqs2[i[1]].id for i in pairs]
    df['seq1_len'] = [len(seqs1[i[0]].seq) for i in pairs]
    df['seq2_len'] = [len(seqs2[i[1]].seq) for i in pairs]
    df['match_len'] = match_lengths
    
    if outfile is not None:
        df.to_csv(outfile,sep='\t')
    
    return df