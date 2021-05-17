# -*- coding: utf-8 -*-
"""
Created on Mon May 17 15:20:21 2021

@author: 20203520
"""


infile_GD = open('GenDescription2.txt')
inlines_GD = infile_GD.read().split()
infile_GD.close()
    
# Voorbeeld_clusterresult.txt inlezen
infile_VC = open('Voorbeeld_clusterresult.txt')
inlines_VC = infile_VC.read().splitlines()
infile_VC.close()



# TWEE DICTIONARIES MAKEN: CLONEID-BESCHRIJVING & CLONEID-CLUSTER
# dictionary CloneID - Beschrijving


indices = [i for i, x in enumerate(inlines_GD) if x == "\\\\"]

beschrijvingdict = {}

for i in range(1,len(indices)):
    line = inlines_GD[indices[i-1]+1: indices[i]]
    ID = int(line[0])
    
    # if i == 0:
    #     Positie[:indices[i]]
    
    Data = ""
    for word in line[1:]:
        if i == 0:
            Positie[:indices[i]]
        Data += " " + word
   
    beschrijvingdict[ID] = Data
    
#print(beschrijvingdict)
# i = 1
# beschrijvingdict = {}

# for i in range(len(inlines_GD)): 
#     lijnen = str(i)
#     lijnen=lijnen.replace('\\\\','')
#     i=int(lijnen)
#     beschrijving = inlines_GD[i]                # beschrijving is regel 2 van het bestand
#     beschrijving = beschrijving.lower()         # alle woorden lowercase
#     beschrijving = beschrijving.split()         # alle woorden in de beschrijving los
#     cloneID = inlines_GD[i-1]                   # van een integer terug naar een string
#     beschrijvingdict[cloneID] = beschrijving    # beschrijving wordt value van key = cloneID
#     i=i+3                                       # 'verpakt' per 3 regels
    
# dictionary CloneID - Cluster
cloneIDdict = {} 
cloneID = []
for line in inlines_VC: 
                   # voor iedere regel in regels van Voorbeeld_clusterresult
    regel = line.split()                        # alle ID's en clusternummer los van elkaar
    cloneID.append(regel[0])                    # voeg het ID van VC (positie 0) samen met het ID van GD
    cloneIDdict[regel[0]]= int(regel[1])        # clusternummer als integer wordt value van de key cloneID



print(beschrijvingdict, cloneIDdict)