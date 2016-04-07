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
#ls_c2o=file_EditedBlast6_to_best_match(f_blast6=f_blast6,f_query=None,f_subject=None,outputfile="C2O.txt",min_identity=95)
