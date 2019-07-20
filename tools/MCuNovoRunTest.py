# -*- coding: utf-8 -*-
"""
Created on Wed Jan 25 16:20:12 2017

@author: k
"""

import MCuNovoGeneSelectorMain
import os
MCuNovofolder = os.path.dirname(os.path.realpath(__file__))
testFolder = MCuNovofolder+'/test/'

f_makerPep=testFolder + 'maker.txt'
f_cuffpep = testFolder + 'cufflinks.txt'
f_denovo1pep =testFolder + 'trinity.txt'
f_denovo2pep= testFolder + 'bridger.txt'

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

f_blast_m2u = testFolder + "O2A.txt"
f_blast_c2u = testFolder + "C2A.txt"
f_blast_d1_2u = testFolder + "T2A.txt"
f_blast_d2_2u = testFolder + "B2A.txt"


f_uniprotpep = testFolder + "reference.txt"
uniqueNoneWithin_c = False,
uniqueNoneWithin_d1 = False,
uniqueNoneWithin_d2 = False,
coverage=0.7
min_length=200
f_protein2gene = testFolder +"peptides2gene.txt"
forceNewName = True
f_outfolder = testFolder

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