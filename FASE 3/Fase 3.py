# -*- coding: utf-8 -*-
"""
Created on Mon May 17 13:06:34 2021

@author: 20202139
"""

# BESTANDEN INLEZEN 
# GenDescription2.txt inlezen
infile_GD = open('GenDescription2.txt')
inlines_GD = infile_GD.read().replace("\\", "").splitlines()
infile_GD.close()
    
# Voorbeeld_clusterresult.txt inlezen
infile_VC = open('Voorbeeld_clusterresult.txt')
inlines_VC = infile_VC.readlines()
infile_VC.close()

    
# TWEE DICTIONARIES MAKEN: CLONEID-BESCHRIJVING & CLONEID-CLUSTER
# dictionary CloneID - Beschrijving
i = 1
beschrijvingdict = {}

for i in range(len(inlines_GD)):                # i moet een regel zijn in de regels van het bestand
    beschrijving = inlines_GD[i]                # beschrijving is regel 2 van het bestand
    beschrijving = beschrijving.split()         # alle woorden in de beschrijving los
    cloneID = str(inlines_GD[i-1])              # regel 0
    beschrijvingdict[cloneID] = beschrijving    # beschrijving wordt value van key = cloneID
    i=i+3                                       # 'verpakt' per 3 regels
    
# dictionary CloneID - Cluster
cloneIDdict = {} 
cloneID = []
for line in inlines_VC:                         # voor iedere regel in regels van Voorbeeld_clusterresult
    regel = line.split()                        # alle ID's en clusternummer los van elkaar
    cloneID.append(regel[0])                    # voeg het ID van VC (positie 0) samen met het ID van GD
    cloneIDdict[regel[0]]= int(regel[1])        # clusternummer als integer wordt value van de key cloneID


print(beschrijvingdict)