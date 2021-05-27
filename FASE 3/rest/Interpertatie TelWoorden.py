# -*- coding: utf-8 -*-
"""
Created on Tue May  4 08:41:37 2021

@author: 20203520
"""


# def Bestand_inlezen(bestandsnaam):
#     infile = open(bestandsnaam)
#     inlines = infile.readlines()
#     infile.close()

# Bestand_inlezen('Voorbeeld_clusterresult.txt')
# Bestand_inlezen('GenDescription2.txt')

# def TelWoorden(bestandsnaam):
#     infile = open(bestandsnaam)
#     inlines = infile.read()
#     infile.close()
    
#     #print(inlines)
#     dic={}
#     inlines.lower()
#     for i in inlines.split()[1::3]:
#         if i in dic:
#             dic[i]= dic[i]+1
#         else:
#             dic[i]=1
#     print(dic)
     
def TelWoorden(bestandsnaam):
    infile = open(bestandsnaam)
    inlines = infile.readlines()
    infile.close()
    
    #print(inlines)
    dic={}
    for line in inlines[1::3]:
        line=line.lower()
        line=line.replace('(','')
        line=line.replace(')','')
        words=line.split()
        for i in words:
            if i in dic:
                dic[i]= dic[i]+1
            else:
                dic[i]=1
    print(dic)         
    print(len(dic.keys()))

TelWoorden('GenDescription2.txt')


def descriptionCloneIDs(bestand1, bestand2):
    
    #open gen_description en voorbeeld_cluster_result bestanden
    inFile = open(bestand1)
    raw_data = inFile.readlines()
    inFile.close
    
    voorbeeld_cluster_result = open(bestand2)
    raw_cluster_data = voorbeeld_cluster_result.readlines()
    voorbeeld_cluster_result.close()
    
    #maak lijst met de beschrijving en overeenkomende clone IDs
    description = [line.strip()for line in raw_data [1::3]]
    cloneIDs = [line.strip() for line in raw_data [0::3]]
   
    #maak lege dictionairy genaamd desciptionCloneIDs
    descriptionCloneIDs = {}
    
    #alle cloneIDs overeenkomend met die beschrijving 
    #all cloneIDs corresponding to that description
    for item in range(len(description)):
        if description[item] not in descriptionCloneIDs:
            descriptionCloneIDs[description[item]] = [cloneIDs[item]]
        else:
            descriptionCloneIDs[description[item]].append[cloneIDs[item]]
   
    #maak lege dictionairy genaamd cluster
    clusterFreq = {}
    
    cluster_data = [items for sublist in [line.strip().split for line in raw_cluster_data] for item in sublist]
    
    cloneIDValues = list(descriptionCloneIDs.values())
    cloneIDKeys = list(descriptionCloneIDs.keys())
    
    for description in range(len(cloneIDsValues[description])):
        if cloneIDValues[description][cloneID] in cluster_data:
            cluster_index = cluster_data.index(cloneIDValues[description][cloneID]) + 1
            cluster_value = int(cluster_data[cluster_index])
            clusterFreq[cloneIDKeys[description]][cluster_value] += 1
            
        else: 
            pass
    print(clusterFreq)
    return clusterFreq

descriptionCloneIDs(bestand1='GenDescription2.txt',bestand2='Voorbeeld_clusterresult.txt')
