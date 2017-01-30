# -*- coding: utf-8 -*-
"""
Created on Wed Jan 25 16:20:12 2017

@author: k
"""

import os
print(os.path.realpath(__file__))
print('Yes')
folder = "F:\\Insects\\Anopheles_gambiae\\"
folderA = folder +"Assemblies\\"
folderC = "F:\\Insects\\Anopheles_gambiae\\MCOT2\\Blast\\"
f_makerPep=folder + "Anopheles-gambiae-PEST_PEPTIDES_AgamP4.3.fa"
f_cuffpep = folderA + "AgCufflinks20160310.fa.transdecoder.pepnonredun"
f_denovo1pep = "F:\\Insects\\Anopheles_gambiae\\TrinityNew\\TrinityConbined\\20160406TrinityCombined98.fa.transdecoder.pepNoWithin"
f_denovo2pep= folderA + "20160314BridgerAll.faUN"
f_blast_m2c = folderC + "O2C.txt"
f_blast_c2m = folderC + "C2O.txt"
f_blast_m2d1 = folderC + "O2T.txt"
f_blast_d1_2m = folderC + "T2O.txt"
f_blast_m2d2 = folderC + "O2B.txt"
f_blast_d2_2m = folderC + "B2O.txt"
f_blast_c2d1 = folderC + "C2T.txt"
f_blast_d1_2c = folderC + "T2C.txt"
f_blast_c2d2 = folderC + "C2B.txt"
f_blast_d2_2c = folderC + "B2C.txt"
f_blast_m2u = folderC + "O2A.txt"
f_blast_c2u = folderC + "C2A.txt"
f_blast_d1_2u = folderC + "T2A.txt"
f_blast_d2_2u = folderC + "B2A.txt"
f_uniprotpep = "F:\\Insects\\uniprotkb_arthropoda"
uniqueNoneWithin_c=True,
uniqueNoneWithin_d1=True,
uniqueNoneWithin_d2=True,
coverage=0.7
min_length=200
f_outfolder = "F:\\Insects\\Anopheles_gambiae\\MCOT2\\"
f_protein2gene = f_outfolder +"peptides2gene.txt"