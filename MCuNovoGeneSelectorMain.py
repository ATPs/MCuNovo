# -*- coding: utf-8 -*-
"""
Created on Thu Apr  7 15:04:25 2016

@author: k
work the same as MCOT, make the code more simple
"""

def file_EditedBlast6_to_best_match(f_blast6,f_query=None,f_subject=None,outputfile=None,min_identity=95,returndic = True):
    """
    Given a filename of Blast6, and query_fasta, subject_fasta file, return a list of with the of all query_ids exists in f_blast6.
    The content in blast6 looks like:
    AGAP000002-PA	TR:Q7QEI4_ANOGA	100.00	1013	0	0	1	1013	1	1013	0.0	 2086
    which includes the query_id, subject_id, identity, alignment_length, mismatches, gap_opens, q_start, q_end, \
    s_start, s_end, e_value, bit_score.
    the returned list should looks like
    (query_id, query_length, subject_id, subject_length, matched_length)
    if f_query, f_subject is not provided, return (query_id, None, subject_id,None, matched_length)
    You can provide f_query and f_subject as the fasta file of query and subject, 
    or f_query and f_subject be dictionary using protein id as key and length as value.
    if returndic is True, return a dictionary, use the query_ad as key, the list element as value
    you can also provide f_blast6 as list of output of readline() commard from fileblast6 (it's easy to combine different list)
    """
    if type(f_blast6) == str:
        with open(f_blast6) as f:
            ls_lines = f.readlines()
    elif type(f_blast6) == list:
        ls_lines = f_blast6
    else:
        print("check the f_blast you provided! not a filename, not a list!")
        return None
    ls_blast6=[] #store each element in blast6 file, convert each of them to proper format
    for line in ls_lines:
        element_list = line.split()
        ls_blast6.append([element_list[0],element_list[1],float(element_list[2]),int(element_list[3]),\
        int(element_list[4]),int(element_list[5]),int(element_list[6]),int(element_list[7]),\
        int(element_list[8]),int(element_list[9]),float(element_list[10]),float(element_list[11])])
    #filter elements in ls_blast6 by min_identity
    ls_blast6Filter =[]
    for ele in ls_blast6:
        if ele[2]>min_identity:
            ls_blast6Filter.append(ele)
    #caluculate matched length for each query subject pair
    ##create a dictionary with query_id as key
    from collections import defaultdict
    dc_blast6 = defaultdict(list)
    for ele in ls_blast6Filter:
        dc_blast6[ele[0]].append(ele)
    ##for each query_id, further create sub-dictionary use subject_id as key
    for ele in dc_blast6:
        templs = dc_blast6[ele]
        dc_blast6[ele] = defaultdict(list)
        for ele2 in templs:
            dc_blast6[ele][ele2[1]].append(ele2)
    ##for each element of element of dc_blast6, calculate the matched_length
    ##dc_blast6, first use query_id as key
    ##then use subject as key. The final value is a list of blast6 line with the same query_id and subject_id
    
    #calculate matched length
    ##create a dictionary the same structure as dc_blast6 to store matched_length
    dc_ml={}
    for ele_query in dc_blast6:
        dc_ml[ele_query] = {}
        for ele_subject in dc_blast6[ele_query]:
            dc_ml[ele_query][ele_subject] = 0
            templs = dc_blast6[ele_query][ele_subject]
            if len(templs) == 1:#if there is only one line of match for a query and subject pair
                tempele = templs[0]
                qstart = tempele[6]
                qend = tempele[7]
                sstart = tempele[8]
                send = tempele[9]
                tempml = min(abs(qend-qstart),abs(send-sstart))+1
            else:
                temp_queryset=set()
                temp_subjectset=set()
                for tempele in templs:
                    qstart = tempele[6]
                    qend = tempele[7]
                    sstart = tempele[8]
                    send = tempele[9]
                    temp_queryset.update(range(qstart,qend+1))
                    temp_subjectset.update(range(sstart,send+1))
                tempml = min(len(temp_queryset),len(temp_subjectset))
            dc_ml[ele_query][ele_subject] = tempml
    
    #select best match pair in dc_ml for each query_id
    ls_bestmatch=[]
    for ele_query in dc_ml:
        tempml = max(dc_ml[ele_query].items(),key=lambda item: item[1])
        ls_bestmatch.append([ele_query, None,tempml[0],None,tempml[1]])
    
    from Bio import SeqIO
    if f_query is not None:
        if type(f_query) != dict:
            dc_querylen={}
            with open(f_query) as f:
                for ele in SeqIO.parse(f,"fasta"):
                    seq = ele.seq
                    if seq[-1] == '*':
                        dc_querylen[ele.id] = len(seq)-1
                    else:
                        dc_querylen[ele.id] = len(seq)
        else:
            dc_querylen = f_query
        for ele in ls_bestmatch:
            ele[1] = dc_querylen[ele[0]]
    if f_subject is not None:
        if type(f_subject) != dict:
            dc_subjectlen={}
            with open(f_subject) as f:
                for ele in SeqIO.parse(f,"fasta"):
                    seq = ele.seq
                    if seq[-1] == '*':
                        dc_subjectlen[ele.id] = len(seq)-1
                    else:
                        dc_subjectlen[ele.id] = len(seq)
        else:
            dc_subjectlen = f_subject
        for ele in ls_bestmatch:
            ele[3] = dc_subjectlen[ele[2]]
    
    if outputfile is not None:
        fout = open(outputfile,"w")
        for ele in ls_bestmatch:
            towrite = ""
            for ele2 in ele:
                towrite = towrite +str(ele2)+"\t"
            towrite = towrite[:-1]+"\n"
            fout.write(towrite)
        fout.close()
    if returndic:
        dc_bestmatch={}
        for ele in ls_bestmatch:
            dc_bestmatch[ele[0]] = ele
        return dc_bestmatch
    else:
        return ls_bestmatch

#folderC = "F:\\Insects\\Anopheles_gambiae\\MCOT2\\Blast\\"
#f_blast6 = folderC + "C2O.txt"
#f_query = "F:\\Insects\\Anopheles_gambiae\\Assemblies\\AgCufflinks20160310.fa.transdecoder.pepnonredun"
#f_subject = "F:\\Insects\\Anopheles_gambiae\\Anopheles-gambiae-PEST_PEPTIDES_AgamP4.3.fa"
#ls_c2o=file_EditedBlast6_to_best_match(f_blast6=f_blast6,f_query=f_query,f_subject=f_subject,outputfile="C2O.txt",min_identity=95)


def openfile2lsFasta(filename,fmt = "fasta"):
    """
    given a file name, return a list of fasta in SeqIO format
    """
    from Bio import SeqIO
    return list(SeqIO.parse(open(filename),fmt))


def openfile2dcFasta(filename, fmt = "fasta"):
    """
    given a filename, return a dictionary of fasta in SeqIO format
    """
    from Bio import SeqIO
    return SeqIO.to_dict(SeqIO.parse(open(filename),fmt))


def lsFasta2dcFasta(lsFasta):
    """
    given a list fasta, return a dictionary
    """
    dcFasta = {}
    for ele in lsFasta:
        dcFasta[ele.id] = ele
    return dcFasta

def fastas2diclen(fastas,remove_star = True,filetype = list):
    """
    given list or dictionary of SeqIO fasta file, return a dictionary, with the id as key, length of seq as value.
    if remove_star is True, if the last letter is *, do not count it.
    """
    dclen ={}
    for ele in fastas:
        if filetype == list:
            seq = str(ele.seq)
            key = ele.id
        elif filetype == dict:
            seq = str(fastas[ele].seq)
            key = ele
        else:
            print("unknown type for fastas")
            return None
        seqlen = len(seq)
        if remove_star:
            if seq[-1] == "*":
                seqlen = seqlen -1
        dclen[key] = seqlen
    return dclen
    
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

#errorMatch("ABCDEFG","ACCDEF")
#errorMatch("AB","ACCDGFF")
#errorMatch("BBCDEFG","ACCDEFG",1)
#errorMatch("ACCDEEGH","ABCDACCDEFGH",1)

        
    
def fasta_within_seq_big_withError(myfasta, error_rate = 0.02,kmerlen = 6):
    """
    myfasta is a list of SeqIO elements
    if a sequence is part of the other, with error_rate allowed, then remove this sequence.
    return a list of non-redundant SeqIO fasta
    """
    import time
    time1 = time.time()
    dickmernum = {} #kmer dic, kmer with its seqs
    for dummyi in range(len(myfasta)):
        seq = str(myfasta[dummyi].seq)
        for i in range(len(seq)+1-kmerlen):
            kmernum = seq[i:i+kmerlen]
            if kmernum not in dickmernum:
                dickmernum[kmernum] = set()
            dickmernum[kmernum].add(dummyi)
    print(time.time()-time1)
    
    time1 = time.time() #change values of dickmernum to list
    for kmernum in dickmernum:
        dickmernum[kmernum] = list(dickmernum[kmernum])
    print(time.time()-time1)
    
    time1 = time.time()
    toremove = set()
    from collections import Counter
    for num1 in range(len(myfasta)):
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


    
def MCuNovoGeneSelectorPep(f_makerPep,f_cuffpep, f_denovo1pep, f_denovo2pep, 
                           f_blast_m2c, f_blast_c2m, 
                           f_blast_m2d1, f_blast_d1_2m, f_blast_m2d2, f_blast_d2_2m,
                           f_blast_c2d1, f_blast_d1_2c, f_blast_c2d2, f_blast_d2_2c,
                           f_blast_m2u, f_blast_c2u, f_blast_d1_2u, f_blast_d2_2u,
                           f_uniprotpep, 
                           uniqueNoneWithin_c=True,
                           uniqueNoneWithin_d1=True,
                           uniqueNoneWithin_d2=True,
                           coverage=0.7, min_length=200,
                           f_protein2gene = None, forceNewName = False,
                           f_outfolder=""):
    """
    the main function.
    """
    import time
    time1 = time.time()
    
    dcpep_m = openfile2dcFasta(f_makerPep) #open maker peptide file
    
    if not uniqueNoneWithin_c:
        print("fasta file of cufflinks is not unique")
        lspep_c = openfile2lsFasta(f_cuffpep)
        lspep_cUniqueNoneWithin = fasta_within_seq_big_withError(lspep_c)
        print("finish generating unique peptides from cufflinks peptide file. ")
        dcpep_c = lsFasta2dcFasta(lspep_cUniqueNoneWithin)
    else:
        dcpep_c = openfile2dcFasta(f_cuffpep)
    
    dcpep_denovo1 = openfile2dcFasta(f_denovo1pep)
    dcpep_denovo2 = openfile2dcFasta(f_denovo2pep)
    dcpep_all = {}
    dcpep_all.update(dcpep_m)
    dcpep_all.update(dcpep_c)
    dcpep_all.update(dcpep_denovo1)
    dcpep_all.update(dcpep_denovo2)
    time2 = time.time()
    print("read in all genes in cufflinks makers, denovo1 and denovo2, time used %.0f"%(time2-time1))
    
    
    def badMatch(ql,sl,ml,coverage=0.7, min_length=200):
        """
        determine whether a match is consider as a bad match.
        if ml/(ql*coverage) + ml/min_length <1, bad match, return True. else, return false
        """
        return ml/(ql*coverage) + ml/min_length <1
        
    #compare maker with cufflinks
    dcML_m2c = file_EditedBlast6_to_best_match(f_blast_m2c,f_makerPep,f_cuffpep)
    mcot_makerKeep =[]##store maker ids that have no good match in cufflinks. 
    for ele in dcpep_m:
        if ele not in dcML_m2c:#maker pep that have no match in cufflinks.
            mcot_makerKeep.append(ele)
    for ele in dcML_m2c:
        query_id, query_len, subject_id,subject_len, matched_length = dcML_m2c[ele]
        if badMatch(query_len, subject_len,matched_length):
            mcot_makerKeep.append(query_id)
    time3 = time.time()
    dcpep_genomebased ={}#store cufflinks unique and genes only found is maker
    for ele in mcot_makerKeep:
        dcpep_genomebased[ele] = dcpep_m[ele]
    for ele in dcpep_c:
        dcpep_genomebased[ele] = dcpep_c[ele]
    print("find maker ids do not have good match or have no mathc in cufflinks. time used %.0f"%(time3-time2))
    print("number maker_unique peptides is %d" %(len(mcot_makerKeep)))
    print("number of genome based peptides is %d\n"%(len(dcpep_genomebased)))
    
    
    #generate protein length dictionary
    dclen_all = fastas2diclen(dcpep_all,True,dict)
    dclen_allAjusted = {}
    for ele in dclen_all:
        if "X" in dcpep_all[ele].seq:
            dclen_allAjusted[ele] = dclen_all[ele] * 0.7
        else:
            dclen_allAjusted[ele] = dclen_all[ele]
#    dclen_uniprot = fastas2diclen(openfile2lsFasta(f_uniprotpep))
    time4 = time.time()
    print("finishing generating dictionary of protein length of this species and Uniprot, time used %.0f"%(time4-time3))
    
    
    #connect genome based peptides to denovo peptides.
    ##get matched length for genome based with denovo peptides
    dcML_m2d1 = file_EditedBlast6_to_best_match(f_blast_m2d1,dclen_all,dclen_all)
    dcML_m2d2 = file_EditedBlast6_to_best_match(f_blast_m2d2,dclen_all,dclen_all)
    dcML_c2d1 = file_EditedBlast6_to_best_match(f_blast_c2d1,dclen_all,dclen_all)
    dcML_c2d2 = file_EditedBlast6_to_best_match(f_blast_c2d2,dclen_all,dclen_all)
    
     ##filter, only keep those matches with ml/ql >coverage (0.7)   
    def bestMatchfilter(dcMatches,coveragelimit = coverage):
        """
        for returned dictionary of file_EditedBlast6_to_best_match, if ml/ql >0.7, keep this match. else, discard
        """
        dcKeep ={}
        for ele in dcMatches:
            query_id, query_len, subject_id,subject_len, matched_length = dcMatches[ele]
            if matched_length/query_len > coveragelimit:
                dcKeep[ele] = dcMatches[ele]
        return dcKeep
    dcML_genomebased2d1={}
    dcML_genomebased2d1.update(bestMatchfilter(dcML_m2d1))
    dcML_genomebased2d1.update(bestMatchfilter(dcML_c2d1))
    dcML_genomebased2d2={}
    dcML_genomebased2d2.update(bestMatchfilter(dcML_m2d2))
    dcML_genomebased2d2.update(bestMatchfilter(dcML_c2d2))
    
    ##connect genome based with denovo. list of list, 
    ##with [genomebased_id, denovo1_id, denovo2_id]
    ##if genomebased have no good match in denovo, use denovoid = genomebased_id
    mcot_genomebased=[]
    for ele in dcpep_genomebased:
        ele1 = ele
        if ele in dcML_genomebased2d1:
            ele2 = dcML_genomebased2d1[ele][2]
        else:
            ele2 = ele
        if ele in dcML_genomebased2d2:
            ele3 = dcML_genomebased2d2[ele][2]
        else:
            ele3 = ele
        mcot_genomebased.append([ele1,ele2,ele3])
    print("finishing matching of maker, cufflinks to denovo1 and denovo2")
    time5 = time.time()
    print("time used %.0f"%(time5-time4))
    
    
    lsBlast2uni=open(f_blast_m2u).readlines() + open(f_blast_c2u).readlines()+ \
    open(f_blast_d1_2u).readlines() + open(f_blast_d2_2u).readlines()
    dcML_all2u = file_EditedBlast6_to_best_match(lsBlast2uni,None,None,None,30,True)#get the matched length of each
    del lsBlast2uni
    def getUniprotlenDsDic(f_uniprotpep,dcML_all2u):
        ls = openfile2lsFasta(f_uniprotpep)
        stUni=set()
        dclen_uniprot ={}
        dcDs_uniprot ={}
        for ele in dcML_all2u:
            stUni.add(dcML_all2u[ele][2])
        for ele in ls:
            if ele.id in stUni:
                dclen_uniprot[ele.id] = len(ele.seq)
                dcDs_uniprot[ele.id] = ele.description
        return dclen_uniprot, dcDs_uniprot
    dclen_uniprot, dcDs_uniprot = getUniprotlenDsDic(f_uniprotpep,dcML_all2u)
            
        
    mcot_genomebased2uniMatchScore=[]
    #for mcot_genomebased, each element contains three ids.
    #get the matching score with uniprot fore each of these ids, 
    #matching score for [genomebased_id, denovo1_id, denovo2_id]
    for ele1 in mcot_genomebased:
        _scores=[]
        for ele2 in ele1:
            if ele2 not in dcML_all2u:
                _scores.append(0)
            else:
                query_id, query_len, subject_id,subject_len, matched_length = dcML_all2u[ele2]
                query_len = dclen_all[query_id]
                subject_len = dclen_uniprot[subject_id]
                _ml = matched_length*matched_length/(query_len*subject_len)
                _scores.append(_ml)
        mcot_genomebased2uniMatchScore.append(_scores)
    time6 = time.time()
    print("finishing calculate match score with uniprot. time used %.0f"%(time6-time5))
    
    
    #select best model in mcot_genomebased
    mcot_genomebasedSelection =[]
    for num in range(len(mcot_genomebased)):
        ele = mcot_genomebased[num] #three ids to select
        elel = [dclen_allAjusted[i] for i in ele] #Ajusted lengthes of these three ids
        elel[0] = elel[0] +5 #to ensure the preference for genome based model, add 5 to genome based length
        eleMatchU = mcot_genomebased2uniMatchScore[num] #match score with uniprot of 3 ids
        if max(eleMatchU) <0.5: # if the best match score with uniprot is less than 0.5. do not use this info, set all matching score to 0
            eleMatchU = [0,0,0]
        elerank =[[ele[i],elel[i],eleMatchU[i]] for i in range(3)]
        #elerank, list of [id, length, matchingscore]
        elerank = sorted(elerank,key= lambda x:x[1])
        #sort by length first. the largest length is the last one of the list
        eleScore = [[elerank[i][0],i + elerank[i][2]*10/3] for i in range(3)] 
        #matching score times 10/3. if the difference is greater than 0.3, it may influce the final rank
        eleScoreRank = sorted(eleScore,key = lambda x:x[1],reverse = True)
        #sort by score. score equals rank + matching score * 10/3. The first element is the best one.
        selection = [ele[0],eleScoreRank[0][0]]
        #print(ele,elel,eleMatchU,selection)
        mcot_genomebasedSelection.append(selection)#[original genome based id, selected id]
    time7 = time.time()
    print("mcot genome based select the best finished. time used %.0f seconds"%(time7-time6))
    
    #find denovo unique proteins.
    dcpep_denovo ={}
    dcpep_denovo.update(dcpep_denovo1)
    dcpep_denovo.update(dcpep_denovo2)
    if (not uniqueNoneWithin_d1) or (not uniqueNoneWithin_d2):
        print("Not using Unique peptides in denovo1 or denovo2!")
        lspep_denovo = list(dcpep_denovo.values())
        lspep_denovoUN = fasta_within_seq_big_withError(lspep_denovo)
        dcpep_denovo ={}
        for ele in lspep_denovoUN:
            dcpep_denovo[ele.id] = ele
    dcML_d1_2c = file_EditedBlast6_to_best_match(f_blast_d1_2c,dclen_all,dclen_all)
    dcML_d1_2m = file_EditedBlast6_to_best_match(f_blast_d1_2m,dclen_all,dclen_all)
    dcML_d2_2c = file_EditedBlast6_to_best_match(f_blast_d2_2c,dclen_all,dclen_all)
    dcML_d2_2m = file_EditedBlast6_to_best_match(f_blast_d2_2m,dclen_all,dclen_all)
    dcML_d2c = {}#denovo to cufflinks
    dcML_d2m = {}#denovo to maker
    dcML_d2c.update(dcML_d1_2c)
    dcML_d2m.update(dcML_d1_2m)
    dcML_d2c.update(dcML_d2_2c)
    dcML_d2m.update(dcML_d2_2m)
    
    dcML_c2m = file_EditedBlast6_to_best_match(f_blast_c2m,dclen_all,dclen_all)
    
    
    
    mcot_denovo =[]##store denovo SeqIO fastas that have no good match in cufflinks and makers. 
    for ele in dcpep_denovo:
        if ele in dcML_d2c:
            query_id, query_len, subject_id,subject_len, matched_length = dcML_d2c[ele]
            denovobadmatch_c = badMatch(query_len, subject_len,matched_length)
        else:
            denovobadmatch_c = True
        if ele in dcML_d2m:
            query_id, query_len, subject_id,subject_len, matched_length = dcML_d2m[ele]
            denovobadmatch_m = badMatch(query_len, subject_len,matched_length)
        else:
            denovobadmatch_m = True
        if denovobadmatch_c and denovobadmatch_m:
            mcot_denovo.append(dcpep_denovo[ele])
    mcot_denovofilter =[]
    mcot_denovobased=[]
    for ele in mcot_denovo:
        if ele.id in dcML_all2u:
            query_id, query_len, subject_id,subject_len, matched_length = dcML_all2u[ele.id]
            query_len = dclen_all[query_id]
            subject_len = dclen_uniprot[subject_id]
            _ml = matched_length*matched_length/(query_len*subject_len) #_ml, matching score
            if _ml>0.6 and "type:complete" in ele.description :
                mcot_denovofilter.append(ele)
            elif "type:complete" in ele.description and len(ele.seq) >200:
                mcot_denovofilter.append(ele)
    mcot_denovoUnique = fasta_within_seq_big_withError(mcot_denovofilter) #remove redundant ones
    for ele in mcot_denovoUnique:
        mcot_denovobased.append([ele.id,ele.id])
    time8 = time.time()
    print("find denovo ids do not have good match or have no match in maker and cufflinks. time used %.0f"%(time8-time7))
    print("number denovo peptides is %d" %(len(mcot_denovo)))
    print("number of filtered denovo peptides is %d\n"%(len(mcot_denovofilter)))
    print("final kept unique denovo peptides is %d"%len(mcot_denovobased))
    
    #finally, get all peptides and annotation to output files.
    mcotSelected = mcot_genomebasedSelection +mcot_denovobased
    mcotSelectedDic={}#use final selected id as key, value [ori,orilen,final, final_len]
    for ele in mcotSelected:
        ele_ori = ele[0]
        ele_final = ele[1]
        ele_orilen = dclen_all[ele_ori]
        ele_finallen = dclen_all[ele_final]
        mcotSelectedDic[ele_final] = [ele_ori,ele_orilen,ele_final,ele_finallen]
    
    mcotpeps=[]
    for ele in mcotSelected:
        mcotpeps.append(dcpep_all[ele[1]])
    print("number of peptides in mcot raw is %d"%len(mcotpeps))
    mcotpepsUN = fasta_within_seq_big_withError(mcotpeps)
    print("number of unique peptides in mcot raw is %d"%len(mcotpepsUN))
    time9 = time.time()
    print("time used %.0f seconds"%(time9-time8))
    
    def return_des(ele,dcMLtemp, dclentemp):
        if ele in dcMLtemp:
            return dcMLtemp[ele][2],dclentemp[dcMLtemp[ele][2]],dcMLtemp[ele][4]
        else:
            return "NA", 0, 0
    
    mcotpepsUNinfo =[]
    for ele in mcotpepsUN:
        orifinal = mcotSelectedDic[ele.id]
        ele_final = orifinal[2]
        ele_finallen = orifinal[3]
        if ele.id in dcpep_m:
            ele_maker, ele_makerlen, ele_makerML = ele.id, ele_finallen, ele_finallen
            ele_cuff,ele_cufflen,ele_cuffML = return_des(ele.id, dcML_m2c,dclen_all)
        elif ele.id in dcpep_c:
            ele_maker, ele_makerlen, ele_makerML = return_des(ele.id, dcML_c2m, dclen_all)
            ele_cuff,ele_cufflen,ele_cuffML = ele.id, ele_finallen, ele_finallen
        if ele.id in dcpep_c or ele.id in dcpep_m:
            ele_denovo1, ele_denovo1len, ele_denovo1ML = return_des(ele.id, dcML_genomebased2d1,dclen_all)
            ele_denovo2, ele_denovo2len, ele_denovo2ML = return_des(ele.id, dcML_genomebased2d2,dclen_all)
        else:
            if ele.id in dcpep_denovo1:
                ele_maker, ele_makerlen, ele_makerML = return_des(ele.id, dcML_d1_2m, dclen_all)
                ele_cuff, ele_cufflen, ele_cuffML = return_des(ele.id, dcML_d1_2c, dclen_all)
                ele_denovo1, ele_denovo1len, ele_denovo1ML  = ele.id, ele_finallen, ele_finallen
                ele_denovo2, ele_denovo2len, ele_denovo2ML = "NA", 0, 0
            else:
                ele_maker, ele_makerlen, ele_makerML = return_des(ele.id, dcML_d2_2m, dclen_all)
                ele_cuff, ele_cufflen, ele_cuffML = return_des(ele.id, dcML_d2_2c, dclen_all)
                ele_denovo1, ele_denovo1len, ele_denovo1ML  = "NA", 0, 0
                ele_denovo2, ele_denovo2len, ele_denovo2ML = ele.id, ele_finallen, ele_finallen
        ele_uni, ele_unilen, ele_uniML = return_des(ele.id, dcML_all2u,dclen_uniprot)
        if ele_uni in dcDs_uniprot:
            ele_uniDes = dcDs_uniprot[ele_uni]
        else:
            ele_uniDes = "NA"
        mcotpepsUNinfo.append(orifinal+[ele_maker, ele_makerlen, ele_makerML]+\
        [ele_cuff,ele_cufflen,ele_cuffML] +[ele_denovo1, ele_denovo1len, ele_denovo1ML]+\
        [ele_denovo2, ele_denovo2len, ele_denovo2ML]+\
        [ele_uni, ele_unilen, ele_uniML,ele_uniDes])
    
    ##assign gene numbers to each of these transcripts.
    mcotpepsUNnames=[]
    if f_protein2gene is None:
        for num in len(mcotpepsUNinfo):
            mcotpepsUNnames[num].append("MCD%07d"%(num+1))
    else:
        dicPep2Gene ={}#store pep2gene information from file
        genelist = open(f_protein2gene).readlines()
        for ele in genelist:
            if '\t' in ele:
                dicPep2Gene[ele.split()[0]] = ele.split()[1]
#        if not forceNewName:
#            dicGene2name ={}
#            for ele in dicPep2Gene:
#                geneid = dicPep2Gene[ele]
#                dicGene2name[geneid] =0
#            genenumbers = len(dicGene2name)
#            for num in range(len(mcotpepsUNinfo)):
#                ele_ori = mcotpepsUNinfo[num][0]
#                if ele_ori in dicPep2Gene:
#                    geneid = dicPep2Gene[ele_ori]
#                    if geneid in dicGene2name:
#                        dicGene2name[geneid] += 1
#                    geneid =geneid+"."+str(dicGene2name[geneid])
#                    mcotpepsUNnames.append(geneid)
#                else:
#                    mcotpepsUNnames.append("MCD%07d"%(genenumbers+1))
#                    genenumbers+=1
#        else:
        dicGene2pepInMcot = {}#geneid as key, list of mcot_ori as element
        for num in range(len(mcotpepsUNinfo)):
            ele_ori = mcotpepsUNinfo[num][0]
            if ele_ori in dicPep2Gene:
                geneid = dicPep2Gene[ele_ori]
                if geneid not in dicGene2pepInMcot:
                    dicGene2pepInMcot[geneid] = []
                dicGene2pepInMcot[geneid].append(ele_ori)
        dicMCOTpep2name = {}#gemerate name for 
        for n, ele in enumerate(dicGene2pepInMcot.items()):
            elekey, elevalue = ele
            for m,c in enumerate(elevalue):
                if forceNewName:
                    if len(elevalue) >1:
                        dicMCOTpep2name[c] = "MCD%07d.%d"%(n+1,m+1)
                    else:
                        dicMCOTpep2name[c] = "MCD%07d.0"%(n+1)
                else:
                    dicMCOTpep2name[c] = elekey+'.'+str(m+1)
        for num in range(len(mcotpepsUNinfo)):
            ele_ori = mcotpepsUNinfo[num][0]
            if ele_ori in dicMCOTpep2name:
                mcotpepsUNnames.append(dicMCOTpep2name[ele_ori])
            else:
                n += 1
                mcotpepsUNnames.append("MCD%07d.0"%(n+1))
            
            
    #save the file of peptides.
    fout_mcot = open(f_outfolder+"MCDpeptides.txt","w")
    to_write =[]
    for num in range(len(mcotpepsUNinfo)):
        to_write.append([">" + mcotpepsUNnames[num] +"\tLen:"+str(mcotpepsUNinfo[num][3])+"\t"+"\t".join([str(i) for i in mcotpepsUNinfo[num]])+"\n",\
        str(mcotpepsUN[num].seq)+"\n"])
    to_write.sort(key=lambda x: x[0])
    for ele in to_write:
        fout_mcot.write(ele[0]+ele[1])
    fout_mcot.close()
    time10 = time.time()
    print("MCOT finished. Total time used %.0f"%(time10-time1))
    print("MCuNovo output:")
    print("MCD gene name; MCD length; oriID, oriLen;finalID,finalLen;")
    print("Maker(ID,len,ML);Cuff(ID,len,ML);denovo1(ID,len,ML);denovo2(ID,len,ML);uniprot(ID,len,ML,description)")
            
            
    
    ##gene number determine function.
    


#folder = "F:\\Insects\\Anopheles_gambiae\\"
#folderA = folder +"Assemblies\\"
#folderC = "F:\\Insects\\Anopheles_gambiae\\MCOT2\\Blast\\"
#f_makerPep=folder + "Anopheles-gambiae-PEST_PEPTIDES_AgamP4.3.fa"
#f_cuffpep = folderA + "AgCufflinks20160310.fa.transdecoder.pepnonredun"
#f_denovo1pep = "F:\\Insects\\Anopheles_gambiae\\TrinityNew\\TrinityConbined\\20160406TrinityCombined98.fa.transdecoder.pepNoWithin"
#f_denovo2pep= folderA + "20160314BridgerAll.faUN"
#f_blast_m2c = folderC + "O2C.txt"
#f_blast_c2m = folderC + "C2O.txt"
#f_blast_m2d1 = folderC + "O2T.txt"
#f_blast_d1_2m = folderC + "T2O.txt"
#f_blast_m2d2 = folderC + "O2B.txt"
#f_blast_d2_2m = folderC + "B2O.txt"
#f_blast_c2d1 = folderC + "C2T.txt"
#f_blast_d1_2c = folderC + "T2C.txt"
#f_blast_c2d2 = folderC + "C2B.txt"
#f_blast_d2_2c = folderC + "B2C.txt"
#f_blast_m2u = folderC + "O2A.txt"
#f_blast_c2u = folderC + "C2A.txt"
#f_blast_d1_2u = folderC + "T2A.txt"
#f_blast_d2_2u = folderC + "B2A.txt"
#f_uniprotpep = "F:\\Insects\\uniprotkb_arthropoda"
#uniqueNoneWithin_c=True,
#uniqueNoneWithin_d1=True,
#uniqueNoneWithin_d2=True,
#coverage=0.7
#min_length=200
#f_outfolder = "F:\\Insects\\Anopheles_gambiae\\MCOT2\\"
#f_protein2gene = f_outfolder +"peptides2gene.txt"

