# -*- coding: utf-8 -*-
"""
Created on Wed May 12 13:05:12 2021

@author: 20202876
"""


def descriptionCloneIDs(gen_description, voorbeeld_cluster_result):
    
    #voor we beginnen, wat willen we clusteren welke data?
    #open gen_description 
    inFile = open(gen_description)
    raw_data = inFile.readlines()
    inFile.close
    
    #open voorbeeld_cluster_result bestanden
    inFile = open(voorbeeld_cluster_result)
    raw_cluster_data = inFile.readlines()
    inFile.close 
    
    #maak lijst met de beschrijving en overeenkomende clone IDs
    description = [[line.strip()] for line in raw_data[1::3]]
    cloneIDs= [[line.strip()] for line in raw_data[0::3]]
   
    #maak lege dictionairy genaamd desciptionCloneIDs
    descriptionCloneIDs = {}
    
    #alle cloneIDs overeenkomend met die beschrijving 
    #doe dit voor het hele bestand
    # for item in range(len(description)):
    #     if description[item] not in descriptionCloneIDs:
    #         descriptionCloneIDs[description[item]] = [cloneIDs[item]]
    #     else:
    #         descriptionCloneIDs[description[item]].append[cloneIDs[item]]
   
    #maak lege dictionairy genaamd cluster
    clusterFreq = {}
    
   # cluster_data = [items for sublist in [line.strip().split and for line in raw_cluster_data] for item in sublist]
    
    cloneIDValues = list(descriptionCloneIDs.values())
    cloneIDKeys = list(descriptionCloneIDs.keys())
    
    # result = []
    # for item in myList:
    #     items = [sub_item.strip() for sub_item in item.split(',')]
    #     items = [ int(sub_item) if sub_item.isdigit() else sub_item for sub_item in items]
    # result.append(items)

    print(len(cloneIDs))       

descriptionCloneIDs('gen_description.txt','voorbeeld_cluster_result.txt')