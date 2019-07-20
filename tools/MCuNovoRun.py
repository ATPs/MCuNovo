# -*- coding: utf-8 -*-
"""
Created on Wed Jan 25 16:20:12 2017

@author: k
"""

import MCuNovoGeneSelectorMain
import os
MCuNovofolder = os.path.dirname(os.path.realpath(__file__))
testFolder = MCuNovofolder+'/test/'



#maker peptide
f_makerPep=testFolder + 'maker.txt'

#cufflinks peptide
f_cuffpep = testFolder + 'cufflinks.txt'

#peptide from de novo program 1
f_denovo1pep =testFolder + 'trinity.txt'

#peptide from de novo program 2
f_denovo2pep= testFolder + 'bridger.txt'





#comparison between 4 different gene models.
#f_blast_x2x, x is m, c, d1, d2, which stand for maker, cufflinks, denovo1, denovo2
#X2X.txt, X is O,C,T,B, which stand for OGS, Cufflinks, Trinity, and Bridger. They are m, c, d1 and d2 shown above
f_blast_m2c = testFolder + "O2C.txt"
f_blast_c2m = testFolder + "C2O.txt"
f_blast_m2d1 = testFolder + "O2T.txt"
f_blast_d1_2m = testFolder + "T2O.txt"
f_blast_m2d2 = testFolder + "O2B.txt"
f_blast_d2_2m = testFolder + "B2O.txt"
f_blast_c2d1 = testFolder + "C2T.txt"
f_blast_d1_2c = testFolder + "T2C.txt"
f_blast_c2d2 = testFolder + "C2B.txt"
f_blast_d2_2c = testFolder + "B2C.txt"



#comparison between modeled genes and reference sequences, here is uniprot_arthropoda seqs. 
#A stands for arthropoda, the reference sequences.
f_blast_m2u = testFolder + "O2A.txt"
f_blast_c2u = testFolder + "C2A.txt"
f_blast_d1_2u = testFolder + "T2A.txt"
f_blast_d2_2u = testFolder + "B2A.txt"

#full path of the reference sequences
f_uniprotpep = testFolder + "reference.txt"





#whether duplicated sequences have been removed from cufflinks, de novo 1 and de novo 2
#set True if you run the optional step to remove duplicated sequences
uniqueNoneWithin_c = False
uniqueNoneWithin_d1 = False
uniqueNoneWithin_d2 = False





#parameters used to determine whether two sequences should be the same. Don't change if you are unclear about the concept
coverage=0.7
min_length=200





#file describe the protein-gene relationship. If None, peptides will be considered from different genes
#set f_protein2gene = None if you do not have peptides2gene.txt file.
f_protein2gene = testFolder +"peptides2gene.txt"

#valid if f_protein2gene is not None. if forceNewName = False, use the gene name provided in the "peptides2gene.txt"
forceNewName = True



#which folder to store the output "MCDpeptides.txt" file
f_outfolder = testFolder





#don't change the code below. Run the main function.
MCuNovoGeneSelectorMain.MCuNovoGeneSelectorPep(f_makerPep,f_cuffpep, f_denovo1pep, f_denovo2pep, 
                           f_blast_m2c, f_blast_c2m, 
                           f_blast_m2d1, f_blast_d1_2m, f_blast_m2d2, f_blast_d2_2m,
                           f_blast_c2d1, f_blast_d1_2c, f_blast_c2d2, f_blast_d2_2c,
                           f_blast_m2u, f_blast_c2u, f_blast_d1_2u, f_blast_d2_2u,
                           f_uniprotpep, 
                           uniqueNoneWithin_c,
                           uniqueNoneWithin_d1,
                           uniqueNoneWithin_d2,
                           coverage, min_length,
                           f_protein2gene, forceNewName,
                           f_outfolder)