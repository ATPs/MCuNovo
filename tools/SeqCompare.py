from multiprocessing import Pool
import sys
sys.path.append(__file__)
import pandas as pd
import time

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
    time1 = time.time()
    print('number of sequences in query', len(seqs1))
    print('number of sequences in subject', len(seqs2))
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
        df.to_csv(outfile,sep='\t', index=None)
    print('done! total seconds: {:.0f}'.format(time.time() - time1))
    return df

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
    parser.add_argument('-q','--query', help = 'a list of SeqIO or a filename of fasta sequence', required=True)
    parser.add_argument('-s','--subject',help = 'a list of SeqIO or a filename of fasta sequence', required=True)
    parser.add_argument('-m','--min_identity',help='the minimum length of identical AAs required to compare two sequences, default 20', type=int, default=20)
    parser.add_argument('-o', '--outfile', help='whether to put output file. default None, only return a dataframe', default=None)
    parser.add_argument('-e', '--error_rate', help='allowed error rate when counts identical AAs, default 0.02', default=0.02, type=float)
    parser.add_argument('-p', '--muscle', help='location of muscle, default = muscle', default='muscle')
    parser.add_argument('-t', '--threads', help='number of CPUs to use, default=8', default=8, type=int)

    f = parser.parse_args()
    compare2lsSeqs(seqs1=f.query, seqs2=f.subject, identitymin = f.min_identity, outfile = f.outfile, error_rate = f.error_rate, muscle_exe = f.muscle, threads=f.threads)
