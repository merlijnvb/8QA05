# -*- coding: utf-8 -*-
"""
Created on Mon May 17 13:06:34 2021

@author: 20202139
"""

# BESTANDEN INLEZEN 
# GenDescription2.txt inlezen
infile_GD = open('GenDescription2.txt')
inlines_GD = infile_GD.read()
infile_GD.close()
    
# Voorbeeld_clusterresult.txt inlezen
infile_VC = open('Voorbeeld_clusterresult.txt')
inlines_VC = infile_VC.read()
infile_VC.close()


# TWEE DICTIONARIES MAKEN: CLONEID-BESCHRIJVING & CLONEID-CLUSTER
# dictionary CloneID - Beschrijving
beschrijvingdict = {}
i = 1
for i in range(len(inlines_GD)):
    beschrijving = inlines_GD[i].split()        # haalt woorden los
    cloneID = inlines_GD[i-1]                # regel 1-1= 0 = de eerste regel
    beschrijvingdict[cloneID] = beschrijving    # geeft de value van cloneID als keys
    i=i+3                                       # per 3 'verpakt' dus dat is de laatste regel
    
# dictionary CloneID - Cluster
cloneID = []
clusterdict = {}   
for line in inlines_VC:
    regel = line.split()                        # haalt de cloneID en clusternummer los (zelfde regel)
    cloneID.append(0)                           # voegt het eerste deel (cloneID) toe
    clusterdict[regel[0]] = int[regel[1]]       # geeft clusternummer als integer en als value




