# -*- coding: utf-8 -*-

import pandas as pd
from Bio import SeqIO
import time
import numpy as np
import math
import fastaFileRemoveDup

pep_m = r"D:\species\genomes\20190622Rhyzopertha_dominica\20190622Rhyzopertha_dominica\Rdo_Scaffolds.all.maker.proteins.fasta"
pep_c = r"D:\species\genomes\20190622Rhyzopertha_dominica\StringTie\20190625RdoStringTie.fa.transdecoder.pepUN"
pep_d = r"D:\species\genomes\20190622Rhyzopertha_dominica\denovo\20190624Trinity.fasta.transdecoder.pepUN"
pep_n = r"D:\species\genomes\20190622Rhyzopertha_dominica\denovo\20190625transabyss.fasta.transdecoder.pepUN"
pep_u = r"D:\species\uniprot\arthropoda\uniprotkb_arthropodaclean"
pair_m2c = r"D:\species\genomes\20190622Rhyzopertha_dominica\20190714MCuNovo\maker2StringTie"
pair_m2d = r"D:\species\genomes\20190622Rhyzopertha_dominica\20190714MCuNovo\maker2Trinity"
pair_m2n = r"D:\species\genomes\20190622Rhyzopertha_dominica\20190714MCuNovo\maker2TransAbyss"
pair_m2u = r"D:\species\genomes\20190622Rhyzopertha_dominica\20190714MCuNovo\Maker2U.match"
pair_c2d = r"D:\species\genomes\20190622Rhyzopertha_dominica\20190714MCuNovo\StringTie2Trinity"
pair_c2n = r"D:\species\genomes\20190622Rhyzopertha_dominica\20190714MCuNovo\StringTie2TransAbyss"
pair_c2u = r"D:\species\genomes\20190622Rhyzopertha_dominica\20190714MCuNovo\StringTie2U.match"
pair_d2u = r"D:\species\genomes\20190622Rhyzopertha_dominica\20190714MCuNovo\Trinity2U.match"
pair_n2u = r"D:\species\genomes\20190622Rhyzopertha_dominica\20190714MCuNovo\TransAbyss2U.match"
protein2transcript2gene = r"D:\species\genomes\20190622Rhyzopertha_dominica\20190714MCuNovo\protein2transcript2gene.txt"
X_table = r"D:\species\genomes\20190622Rhyzopertha_dominica\20190714MCuNovo\proteinXcount.txt"
outprefix = r'D:\species\genomes\20190622Rhyzopertha_dominica\20190714MCuNovo\MCD'
CS_coverage = 0.7
CS_length=200
MS_min=0.5
MS_factor = 10/3
X_ratio = 0.7
homolog_CSmin = 0.6
mc_extra = 5

def openfile2lsFasta(filename,fmt = "fasta", removeStar = True):
    """
    given a file name, return a list of fasta in SeqIO format
    """
    l = list(SeqIO.parse(open(filename),fmt))
    if removeStar:
        for s in l:
            if s.seq[-1] == '*':
                s.seq = s.seq[:-1]
    return l

def isna(x):
    '''
    return True or False, if x is np.nan
    '''
    if isinstance(x,float):
        return np.isnan(x)
    return False

def readPair(filename,sep='\t'):
    '''
    read in pair comparison file and add MS column
    '''
    df = pd.read_csv(filename, sep='\t')
    df['MS'] = df['match_len'] ** 2 /df['seq1_len'] / df['seq2_len']
    return df

def selectBasedOnDf_mc2dnu(r):
    '''
    each row is like
mc                            MSTRG.10.1.p1
mc_len                                  432
d                TRINITY_DN5924_c0_g1_i1.p1
d_len                                   432
mc2d_ML                                 432
mc2d_MS                                   1
n                              R21486350.p1
n_len                                   432
mc2n_ML                                 432
mc2n_MS                                   1
u                       TR:A0A2J7QRR2_9NEOP
u_len                                   427
mc2u_ML                                 427
mc2u_MS                            0.988426
du                      TR:A0A067RDJ7_ZOONE
du_len                                  427
d2u_ML                                  426
d2u_MS                             0.983802
nu                      TR:A0A067RDJ7_ZOONE
nu_len                                  427
n2u_ML                                  426
n2u_MS                             0.983802
mc_Xcount                                 0
mc_len_adjust                           432
mOnly                                 False
    select based on each row
    start compare if n_len / mc_len > CS_coverage
    return annotation and selected ID
    '''
    #check if the query is from cufflinks or maker
    if r['mOnly']:
        q = 'm'
    else:
        q = 'c'
    #reshape r to a dataframe
    tdf = pd.DataFrame()
    tdf['type'] = [q,'d','n']
    tdf['id'] = [r['mc'], r['d'], r['n']]
    tdf['len'] = [r['mc_len_adjust'] + mc_extra, r['d_len'], r['n_len']] #change mc_len
    tdf['u'] = [r['u'],r['du'],r['nu']]
    tdf['u_len'] = [r['u_len'], r['du_len'], r['nu_len']]
    tdf['u_ML'] = [r['mc2u_ML'], r['d2u_ML'], r['n2u_ML']]
    tdf['u_MS'] = [r['mc2u_MS'], r['d2u_MS'], r['n2u_MS']]
    tdf.index = tdf['type']
    #only keep rows with len/mc_len > CS_coverage
    row_keep = [q]
    if r['d_len'] / r['mc_len'] >= CS_coverage:
        row_keep.append('d')
    if r['n_len'] / r['mc_len'] > CS_coverage:
        row_keep.append('n')
    tdf = tdf.loc[row_keep]
    
    #score based on length
    tdf = tdf.sort_values(by='len',ascending=False)
    tdf['scores'] = [3,2,1][:tdf.shape[0]]#assign score based on length
    #add score based on u_MS, if max(u_MS) > homolog_CSmin, then scores = scores + MS_factor * u_MS
    if tdf['u_MS'].max() > homolog_CSmin:
        tdf['scores'] = tdf['scores'] + MS_factor * tdf['u_MS']
    
    #return the selected one
    tdf = tdf.sort_values(by='scores', ascending=False)
    keep = tdf.iloc[0]
    return q, keep['type'], keep['id'], keep['len']

def selector(pep_m, pep_c, pep_d, pep_n, pep_u,
             pair_m2c, pair_m2d, pair_m2n, pair_m2u, pair_c2d, pair_c2n, pair_c2u, pair_d2u, pair_n2u,
             protein2transcript2gene,
             outprefix = '',
             X_table = None,
             CS_coverage = 0.7,
             CS_length=200,
             MS_min=0.5,
             MS_factor = 10/3,
             X_ratio = 0.7,
             homolog_CSmin = 0.6,
             mc_extra = 5,
             error_rate = 0.01
             ):
    '''
    m: maker
    c: cufflinks or StringTie
    d, n: two de novo
    pep_m, pep_c, pep_d, pep_n, pep_u: protein sequences
    pair_m2c, pair_m2d, pair_m2n, pair_m2u, pair_c2d, pair_c2n, pair_c2u, pair_d2u, pair_c2u: one to one best match between different protein sets
    protein2transcript2gene: a table with protein, transcript and gene ids
    coverage: threshold to compare m2c, m/c to d/n
    CS_coverage, threshold to calculate "Confidence Score"
    CS_length: consider "Matching Score" if there is at least one MS greater than CS_length
    MS_factor: factor of contribution of MS when doing selection
    X_table: a table of proteins with X and the Xcount
    X_ratio: how to punish sequences with 'X'. if X_ratio >= 1, adjusted_protein_len = protein_len - XCount * X_ratio; if X_ratio < 1, adjusted_protein_len = protein_len * x_ratio
    outprefix: prefix for output files
    homolog_CSmin: for makerOnly and de novo Only sequences, require them to have homolog with CS >= CS_homolog
    mc_extra: add mc_extra to the length of maker/cufflinks when do selection
    error_rate: allowed error rate when removing duplicated sequences
    '''
    time0 = time.time()
    
    #read in all comparision pairs
    df_m2c = readPair(pair_m2c, sep='\t')
    df_m2d = readPair(pair_m2d, sep='\t')
    df_m2n = readPair(pair_m2n, sep='\t')
    df_m2u = readPair(pair_m2u, sep='\t')
    df_c2d = readPair(pair_c2d, sep='\t')
    df_c2n = readPair(pair_c2n, sep='\t')
    df_c2u = readPair(pair_c2u, sep='\t')
    df_d2u = readPair(pair_d2u, sep='\t')
    df_n2u = readPair(pair_n2u, sep='\t')
    print('time {T:.0f} s: '.format(T=time.time() - time0), end='')
    print('finish readin all comparison pairs')
    
    # compare maker with cufflinks/StringTie
    print('time {T:.0f} s: '.format(T=time.time() - time0), end='')
    print('step1, compare maker with StringTie/Cufflinks'.format(T=time.time() - time0))
    ## read in maker sequences
    ls_seqs_maker = openfile2lsFasta(pep_m)
    print('time {T:.0f} s: finish read in maker sequences. There are {seqcount} proteins in maker'.format(T=time.time() - time0, seqcount = len(ls_seqs_maker)))
    ## calculate Confidence Score
    df_m2c['CS'] = df_m2c['match_len'] / (CS_coverage * df_m2c['seq1_len']) + df_m2c['match_len'] / CS_length
    ## sort by CS and drop duplicates based on seq1_id
    df_m2c = df_m2c.sort_values(by=['seq1_id', 'CS'], ascending=[True, False])
    df_m2c = df_m2c.drop_duplicates(subset='seq1_id', keep='first')
    df_m_keep = df_m2c[df_m2c['CS'] < 1].copy()
    ## keep maker proteins with no match in cufflinks
    st_maker = set(df_m2c['seq1_id'])
    ls_m_unmap = [e for e in ls_seqs_maker if e.id not in st_maker] #store protein id and length
    pep_m_only = set([e.id for e in ls_m_unmap] + list(df_m_keep['seq1_id']))
    ## save sequences
    fout = open(outprefix+'.makerOnly.protein.txt','w')
    for s in ls_seqs_maker:
        if s.id in pep_m_only:
            fout.write('>'+s.description+'\n'+str(s.seq)+'\n')
    fout.close()
    print('time {T:.0f} s: '.format(T=time.time() - time0), end='')
    print('finish read in maker2cufflink comparison. Maker sequences, {count1} with match in cufflinks, of which {count2} with CS < 1. {count3} sequences with no match in cufflinks'.format(count1 = len(st_maker), count2 = df_m_keep.shape[0], count3 = len(ls_m_unmap)))
    ## for makerOnly sequences, compare with uniprot and only keep those with MS > homolog_CSmin
    df_m2u_good = df_m2u[df_m2u['MS'] > homolog_CSmin]
    maker_withGoodU = set(df_m2u_good['seq1_id'])
    pep_m_onlyGood = set([e for e in pep_m_only if e in maker_withGoodU])
    print('of the {count1} maker-only sequences, after filter based on matching with homologs, {count2} maker-only sequences left'.format(count1 = len(pep_m_only),count2 = len(pep_m_onlyGood)))
    ls_m_onlyGood = [e for e in ls_seqs_maker if e.id in pep_m_onlyGood]
    
    # combine StringTie/Cufflinks with makerOnly and compare with de novo
    ## read in cufflinks sequences to list
    ls_seqs_cufflinks = openfile2lsFasta(pep_c)
    ## create df_mc2d, df_mc2n, df_mc2u, each seq_1 is all sequences in StringTies/Cuffinks + pep_m_onlyGood
    ### create df_mc2d
    df_m2d_mOnlyGood = df_m2d[df_m2d['seq1_id'].apply(lambda x:x in pep_m_onlyGood)]
    df_mc2d = pd.concat([df_m2d_mOnlyGood, df_c2d]) #combine df_m2d, df_c2d
    df_mc2d = df_mc2d.sort_values(by=['seq1_id','match_len','seq2_len'], ascending=[True,False,False])
    df_mc2d = df_mc2d.drop_duplicates(subset='seq1_id', keep='first').copy()
    ### add pep_m_onlyGood and ls_seqs_cufflinks with no match in d
    st_mc = set(df_mc2d['seq1_id'])
    ls_mc_noMatchD = [e for e in ls_m_onlyGood if e.id not in st_mc] + [e for e in ls_seqs_cufflinks if e.id not in st_mc]
    df_mc2d_supp = pd.DataFrame(columns=df_mc2d.columns)
    df_mc2d_supp['seq1_id'] = [e.id for e in ls_mc_noMatchD]
    df_mc2d_supp['seq1_len'] = [len(e.seq) for e in ls_mc_noMatchD]
    df_mc2d = pd.concat([df_mc2d, df_mc2d_supp])
    
    ### create df_mc2n
    df_m2n_mOnlyGood = df_m2n[df_m2n['seq1_id'].apply(lambda x:x in pep_m_onlyGood)]
    df_mc2n = pd.concat([df_m2n_mOnlyGood, df_c2n]) #combine df_m2n, df_cnd
    df_mc2n = df_mc2n.sort_values(by=['seq1_id','match_len','seq2_len'], ascending=[True,False,False])
    df_mc2n = df_mc2n.drop_duplicates(subset='seq1_id', keep='first').copy()
    ### add pep_m_onlyGood and ls_seqs_cufflinks with no match in d
    st_mc = set(df_mc2n['seq1_id'])
    ls_mc_noMatchD = [e for e in ls_m_onlyGood if e.id not in st_mc] + [e for e in ls_seqs_cufflinks if e.id not in st_mc]
    df_mc2n_supp = pd.DataFrame(columns=df_mc2n.columns)
    df_mc2n_supp['seq1_id'] = [e.id for e in ls_mc_noMatchD]
    df_mc2n_supp['seq1_len'] = [len(e.seq) for e in ls_mc_noMatchD]
    df_mc2n = pd.concat([df_mc2n, df_mc2n_supp])
    
    ### create df_mc2u
    df_m2u_mOnlyGood = df_m2u[df_m2u['seq1_id'].apply(lambda x:x in pep_m_onlyGood)]
    df_mc2u = pd.concat([df_m2u_mOnlyGood, df_c2u]) #combine df_m2n, df_cnd
    df_mc2u = df_mc2u.sort_values(by=['seq1_id','MS'], ascending=[True,False])
    df_mc2u = df_mc2u.drop_duplicates(subset='seq1_id', keep='first').copy()
    ### add pep_m_onlyGood and ls_seqs_cufflinks with no match in d
    st_mc = set(df_mc2u['seq1_id'])
    ls_mc_noMatchD = [e for e in ls_m_onlyGood if e.id not in st_mc] + [e for e in ls_seqs_cufflinks if e.id not in st_mc]
    df_mc2u_supp = pd.DataFrame(columns=df_mc2u.columns)
    df_mc2u_supp['seq1_id'] = [e.id for e in ls_mc_noMatchD]
    df_mc2u_supp['seq1_len'] = [len(e.seq) for e in ls_mc_noMatchD]
    df_mc2u = pd.concat([df_mc2u, df_mc2u_supp])
    
    ### combine df_mc2d, df_mc2n, df_mc2u
    df_mc2dnu = pd.DataFrame(columns = ['mc', 'mc_len', 'd', 'd_len', 'mc2d_ML', 'mc2d_MS', 'n', 'n_len', 'mc2n_ML', 'mc2n_MS', 'u', 'u_len', 'mc2u_ML', 'mc2u_MS','du', 'du_len', 'd2u_ML', 'd2u_MS', 'nu', 'nu_len', 'n2u_ML', 'n2u_MS'])
    df_mc2dnu['mc'] = df_mc2d['seq1_id']
    df_mc2dnu['mc_len'] = df_mc2d['seq1_len']
    df_mc2dnu['d'] = df_mc2d['seq2_id']
    df_mc2dnu['d_len'] = df_mc2d['seq2_len']
    df_mc2dnu['mc2d_ML'] = df_mc2d['match_len']
    df_mc2dnu['mc2d_MS'] = df_mc2d['MS']
    dc = dict(zip(df_mc2n['seq1_id'], df_mc2n['seq2_id']))
    df_mc2dnu['n'] = df_mc2dnu['mc'].apply(lambda x:dc[x])
    dc = dict(zip(df_mc2n['seq1_id'], df_mc2n['seq2_len']))
    df_mc2dnu['n_len'] = df_mc2dnu['mc'].apply(lambda x:dc[x])
    dc = dict(zip(df_mc2n['seq1_id'], df_mc2n['match_len']))
    df_mc2dnu['mc2n_ML'] = df_mc2dnu['mc'].apply(lambda x:dc[x])
    dc = dict(zip(df_mc2n['seq1_id'], df_mc2n['MS']))
    df_mc2dnu['mc2n_MS'] = df_mc2dnu['mc'].apply(lambda x:dc[x])
    dc = dict(zip(df_mc2u['seq1_id'], df_mc2u['seq2_id']))
    df_mc2dnu['u'] = df_mc2dnu['mc'].apply(lambda x:dc[x])
    dc = dict(zip(df_mc2u['seq1_id'], df_mc2u['seq2_len']))
    df_mc2dnu['u_len'] = df_mc2dnu['mc'].apply(lambda x:dc[x])
    dc = dict(zip(df_mc2u['seq1_id'], df_mc2u['match_len']))
    df_mc2dnu['mc2u_ML'] = df_mc2dnu['mc'].apply(lambda x:dc[x])
    dc = dict(zip(df_mc2u['seq1_id'], df_mc2u['MS']))
    df_mc2dnu['mc2u_MS'] = df_mc2dnu['mc'].apply(lambda x:dc[x])
    ###get bett MS for each d
    df_d2u_best = df_d2u.sort_values(by=['seq1_id','MS'], ascending=[True, False])
    df_d2u_best = df_d2u.drop_duplicates(subset='seq1_id',keep='first')
    dc = dict(zip(df_d2u_best['seq1_id'], df_d2u_best['seq2_id']))
    df_mc2dnu['du'] = df_mc2dnu['d'].apply(lambda x:dc[x] if x in dc else np.nan)
    dc = dict(zip(df_d2u_best['seq1_id'], df_d2u_best['seq2_len']))
    df_mc2dnu['du_len'] = df_mc2dnu['d'].apply(lambda x:dc[x] if x in dc else np.nan)
    dc = dict(zip(df_d2u_best['seq1_id'], df_d2u_best['match_len']))
    df_mc2dnu['d2u_ML'] = df_mc2dnu['d'].apply(lambda x:dc[x] if x in dc else np.nan)
    dc = dict(zip(df_d2u_best['seq1_id'], df_d2u_best['MS']))
    df_mc2dnu['d2u_MS'] = df_mc2dnu['d'].apply(lambda x:dc[x] if x in dc else np.nan)
    ###get bett MS for each n
    df_n2u_best = df_n2u.sort_values(by=['seq1_id','MS'], ascending=[True, False])
    df_n2u_best = df_n2u.drop_duplicates(subset='seq1_id',keep='first')
    dc = dict(zip(df_n2u_best['seq1_id'], df_n2u_best['seq2_id']))
    df_mc2dnu['nu'] = df_mc2dnu['n'].apply(lambda x:dc[x] if x in dc else np.nan)
    dc = dict(zip(df_n2u_best['seq1_id'], df_n2u_best['seq2_len']))
    df_mc2dnu['nu_len'] = df_mc2dnu['n'].apply(lambda x:dc[x] if x in dc else np.nan)
    dc = dict(zip(df_n2u_best['seq1_id'], df_n2u_best['match_len']))
    df_mc2dnu['n2u_ML'] = df_mc2dnu['n'].apply(lambda x:dc[x] if x in dc else np.nan)
    dc = dict(zip(df_n2u_best['seq1_id'], df_n2u_best['MS']))
    df_mc2dnu['n2u_MS'] = df_mc2dnu['n'].apply(lambda x:dc[x] if x in dc else np.nan)
    
    #read in proteinXcount and add 'mc_Xcount' column for df_mc2dnu
    df_Xcount = pd.read_csv(X_table, sep='\t')
    dc_pr2X = dict(zip(df_Xcount['protein'], df_Xcount['Xcount']))
    df_mc2dnu['mc_Xcount'] = df_mc2dnu['mc'].apply(lambda x:dc_pr2X[x] if x in dc_pr2X else 0)
    
    #add 'mc_len_adjust' column for df_mc2dnu
    if X_ratio >=1:
        df_mc2dnu['mc_len_adjust'] = df_mc2dnu['mc_len'] - X_ratio * df_mc2dnu['mc_Xcount']
    elif 0 <= X_ratio and X_ratio <1:
        df_mc2dnu['mc_len_adjust'] = df_mc2dnu.apply(lambda x:x['mc_len'] * X_ratio if x['mc_Xcount'] > 0 else x['mc_len'], axis=1)
    else:
        print('something wrong! for X_ration, cannot be less than 0!')
        exit(1)
    
    # add column 'mOnly' to indicate in the mc column it is a maker id
    df_mc2dnu['mOnly'] = df_mc2dnu['mc'].apply(lambda x:x in pep_m_onlyGood)
    
    # make selection for each row of df_mc2dnu, add Notes
    keeps = df_mc2dnu.apply(selectBasedOnDf_mc2dnu, axis=1)
    
    # add select results to df_mc2dnu
    df_mc2dnu['s_q'] = [e[0] for e in keeps]
    df_mc2dnu['s_s'] = [e[1] for e in keeps]
    df_mc2dnu['s_id'] = [e[2] for e in keeps]
    df_mc2dnu['s_len'] = [e[3] for e in keeps]
    # save df_mc2dnu
    df_mc2dnu.to_csv(outprefix+'.mc2dnu.tab.txt', sep='\t',index=None)
    
    
    # get proteins only in de novo
    ## read in de novo d sequences
    ls_seqs_novo = openfile2lsFasta(pep_d) + openfile2lsFasta(pep_n)
    ## calculate Confidence Score
    df_d2mc = pd.concat([df_m2d, df_m2n, df_c2d, df_c2n])
    df_d2mc.columns = ['seq2_id', 'seq1_id', 'seq2_len', 'seq1_len', 'match_len', 'MS']
    df_d2mc['CS'] = df_d2mc['match_len'] / (CS_coverage * df_d2mc['seq1_len']) + df_d2mc['match_len'] / CS_length
    ## sort by CS and drop duplicates based on seq1_id
    df_d2mc = df_d2mc.sort_values(by=['seq1_id', 'CS'], ascending=[True, False])
    df_d2mc = df_d2mc.drop_duplicates(subset='seq1_id', keep='first')
    df_novo_keep = df_d2mc[df_d2mc['CS'] < 1].copy()#de novo with no good match in maker/cufflinks
    ## keep maker proteins with no match in cufflinks
    st_novo = set(df_d2mc['seq1_id'])
    ls_novo_unmap = [e for e in ls_seqs_novo if e.id not in st_novo] #store protein id and length
    pep_novo_only = set([e.id for e in ls_novo_unmap] + list(df_novo_keep['seq1_id']))
    ## save sequences
    fout = open(outprefix+'.denovoOnly.protein.txt','w')
    for s in ls_seqs_novo:
        if s.id in pep_novo_only:
            fout.write('>'+s.description+'\n'+str(s.seq)+'\n')
    fout.close()
    ## keep de novo with good match in Uniprot
    df_novo2u = pd.concat([df_d2u, df_n2u])
    df_novo2u_good = df_novo2u[df_novo2u['MS'] > homolog_CSmin]
    novo_withGoodU = set(df_novo2u_good['seq1_id'])
    pep_novo_onlyGood = set([e for e in pep_novo_only if e in novo_withGoodU])
    print('of the {count1} denovo-only sequences, after filter based on matching with homologs, {count2} denovo-only sequences left'.format(count1 = len(pep_novo_only),count2 = len(pep_novo_onlyGood)))
    ls_novo_onlyGood = [e for e in ls_seqs_novo if e.id in pep_novo_onlyGood]
    
    ## get df_mc2dnu sequences
    ls_seqs_mc2dnu = []
    ls_seqs_cuff = openfile2lsFasta(pep_c)
    ls_seqs_all = ls_seqs_cuff + ls_seqs_maker + ls_seqs_novo
    mcd_pep_ids_mcdnu = set(list(df_mc2dnu['s_id']))
    for seq in ls_seqs_all:
        if seq.id in mcd_pep_ids_mcdnu:
            ls_seqs_mc2dnu.append(seq)
    
    #remove duplicated sequences
    ls_seqs_mcdall = ls_seqs_mc2dnu + ls_novo_onlyGood
    ls_seqs_mcdUN = fastaFileRemoveDup.fasta_within_seq_big_withError(ls_seqs_mcdall, error_rate=error_rate)
    #filter df_mc2dnu to keep only those in ls_seqs_mcdUN
    mcd_pep_ids = set([e.id for e in ls_seqs_mcdUN])
    df_mc2dnu2 = df_mc2dnu[df_mc2dnu['s_id'].isin(mcd_pep_ids)]
    
    #rename selected sequences
    ## read in protein2transcript2gene
    df_p2g = pd.read_csv(protein2transcript2gene,sep='\t')
    dc_p2g = dict(zip(df_p2g['protein'], df_p2g['gene']))
    ## create df store selected id and their name
    df_name = pd.DataFrame()
    df_name['id'] = list(mcd_pep_ids)
    df_name['gene'] = [str(i) for i in df_name.index]
    df_name['gene'] = df_name.apply(lambda x:dc_p2g[x['id']] if x['id'] in dc_p2g else x['gene'], axis=1)
    df_name['gene_ori'] = df_name['gene']#to compare the number after merging genes
    #merge genes based on matching in df_mc2dnu
    for row,r in df_mc2dnu2.iterrows():
        ids = [r['mc']]
        if r['d_len'] / r['mc_len'] >= CS_coverage:
            ids.append(r['d'])
        if r['n_len'] / r['mc_len'] > CS_coverage:
            ids.append(r['n'])
        ids = [e for e in ids if e in mcd_pep_ids]
        if len(ids) <= 1:
            continue
        else:
            genes = list(df_name[df_name['id'].isin(ids)]['gene'])
            gene = genes[0]
            genes = set(genes)
            if len(genes) <= 1:
                continue
            else:
                df_name.loc[df_name['gene'].isin(genes),'gene'] = gene
                #print(genes)
    print('before merging, there are {} genes; after merging based on matching, there are {} genes'.format(len(df_name['gene_ori'].unique()), len(df_name['gene'].unique())))
    ## assign MCD names
    genes = df_name['gene'].value_counts()
    gene_one = set(genes[genes==1].index) #genes with only one sequence
    df_name = df_name.sort_values(by='gene')
    genes = df_name['gene'].drop_duplicates()
    dc_gene2geneN = dict(zip(genes, range(1,len(genes)+1)))
    df_name['name_gene'] = df_name['gene'].apply(lambda x:dc_gene2geneN[x])
    gene_transcript = []
    genes = list(df_name['gene'])
    if genes[0] in gene_one:
        gene_transcript.append(0)
    else:
        gene_transcript.append(1)
    for i in range(1, len(genes)):
        if genes[i] in gene_one:
            gene_transcript.append(0)
        else:
            if genes[i] == genes[i-1]:
                gene_transcript.append(gene_transcript[-1]+1)
            else:
                gene_transcript.append(1)
    df_name['name_transcript'] = gene_transcript
    df_name['name_mcd'] = df_name.apply(lambda x:'MCD{}.{}'.format(x['name_gene'], x['name_transcript']), axis=1)
    dc_id2name = dict(zip(df_name['id'], df_name['name_mcd']))
    
    #write annotation and peptide sequences
    fout = open(outprefix+'.MCD.protein.withNotes.txt','w')
    fout2 = open(outprefix+'.MCD.protein.txt','w')
    ids_novoOnly = set([e.id for e in ls_novo_onlyGood])
    for seq in ls_seqs_mcdUN:
        if seq.id in ids_novoOnly:#write novo Only sequences
            gene_id = dc_id2name[seq.id]
            fout.write('>{}.novo\n{}\n'.format(gene_id, str(seq.seq)))
            fout2.write('>{}.novo\n{}\n'.format(gene_id, str(seq.seq)))
        else:
            gene_id = dc_id2name[seq.id]
            tdf = df_mc2dnu2[df_mc2dnu2['s_id'] == seq.id]
            if tdf.shape[0] == 1:
                r = tdf.iloc[0]
                gene_id = gene_id + '.{}{}'.format(r['s_q'], r['s_s'])
                description = str(r.to_dict())
                fout.write('>{} {}\n{}\n'.format(gene_id, description, str(seq.seq)))
                fout2.write('>{}\n{}\n'.format(gene_id, str(seq.seq)))
            else:
                description = []
                select_notes = []
                for row, r in tdf.iterrows():
                    description.append(str(r.to_dict()))
                    select_notes.append('{}{}'.format(r['s_q'], r['s_s']))
                select_notes = list(set(select_notes))
                select_notes = '.'.join(select_notes)
                description = '||'.join(description)
                gene_id = gene_id + '.'+select_notes
                fout.write('>{} {}\n{}\n'.format(gene_id, description, str(seq.seq)))
                fout2.write('>{}\n{}\n'.format(gene_id, str(seq.seq)))
    fout.close()
    fout2.close()

    
    




















    
    

