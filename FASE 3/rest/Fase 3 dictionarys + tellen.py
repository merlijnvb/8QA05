# -*- coding: utf-8 -*-
"""
Created on Mon May 17 15:20:21 2021

@author: 20203520
"""
import numpy as np

def lib_res(Filename):
    def openfile(Filename):
        read_data = open(f'{Filename}')
        data = read_data.read()
        data = data.splitlines()
        read_data.close()
        
        return data
    
    lib = dict()
    
    for lines in openfile(Filename):
        line = lines.split()
        lib[int(line[0])] = int(line[1])
        
    return lib

def lib_beschrijvingen(Filename):
    infile_GD = open(f'{Filename}')
    inlines_GD = infile_GD.read().split()
    infile_GD.close() 
    
    indices = [i for i, x in enumerate(inlines_GD) if x == "\\\\"]  # als item is \\\\, bewaar de index
    beschrijvingdict = {}
    
    
    line = inlines_GD[: indices[0]]       # elke lijn is het cloneID met zijn beschrijving
    cloneID = int(line[0])                              # geeft cloneID als integer
    
    beschrijving = ""
    for word in range(1, len(line[1:])+1):   
        if word == 1:
            beschrijving += line[word] 
        
        else:
            beschrijving += " " + line[word]                              # print woorden met spatie
    beschrijvingdict[cloneID] = beschrijving
    
    for i in range(1,len(indices)):   
        line = inlines_GD[indices[i-1]+1: indices[i]]       # elke lijn is het cloneID met zijn beschrijving
        cloneID = int(line[0])                              # geeft cloneID als integer
        
        beschrijving = ""
        for word in range(1, len(line[1:])+1):   
            if word == 1:
                beschrijving += line[word] 

            else:
                beschrijving += " " + line[word]                              # print woorden met spatie
        
        beschrijvingdict[cloneID] = beschrijving
        
    return beschrijvingdict

def get_descrision(input_data, refrence_data):
    lib_refrence = dict()
    
    for ID in input_data:
        lib_refrence[ID] = refrence_data[ID]
        
    return lib_refrence

def get_cluster_descrision(results, refrence_data):
    cluster_discription = dict()
    
    for i in range(min(results.values()),max(results.values())+1):
        cluster_discription[i] = list()
        
    for ID in results:
        cluster_discription[results[ID]].append(refrence_data[ID])
        
    return cluster_discription

def get_same_descriptions(lib_cluster_info):
    clusters = list(lib_cluster_info.keys())
    same_descriptions = {}
    for cluster in clusters:
        descriptions = list(lib_cluster_info[cluster])
        multiple_desc = {}
        for description in descriptions:
            duplicates = any(descriptions.count(description) > 1 for element in descriptions)
            if duplicates == True:
                multiple_desc[description] = descriptions.count(description)
        same_descriptions[cluster] = multiple_desc
        
    return same_descriptions

def count_descriptions(same_descriptions, clusters=1):
    all_clusters = list(same_descriptions.keys())
    
    for cluster in all_clusters:
        cluster_descriptions = list(same_descriptions[cluster].keys())
        alle_clusters = list(same_descriptions.keys())
        alle_clusters.pop(cluster)
        
        for other_cluster in alle_clusters:
            other_cluster_desc = list(same_descriptions[other_cluster].keys())
            
            for desc in cluster_descriptions:
                if any(other_cluster_desc.count(desc) > 0 for element in other_cluster_desc) == True:
                    print(desc)
                    '''nu geeft ie als resultaat alleen ESTs, maar zo als het in de casus staat zou
                    die meer moeten returnen. De fout zit waarschijnlijk in de functie get_same_descriptions,
                    omdat deze nu kijkt naar of de strings precies het zelfde zijn. Maar eigenlijk moet
                    deze functie kijken of de substring in een andere  string van die lijst voorkomt.'''
                
        all_clusters = list(same_descriptions.keys())

#count_descriptions = count_descriptions(same_descriptions)

lib_results = lib_res('Voorbeeld_clusterresult.txt')
lib_beschrijving = lib_beschrijvingen('GenDescription2.txt')
lib_bes_toegekend = get_descrision(lib_results,lib_beschrijving)
lib_cluster_info = get_cluster_descrision(lib_results, lib_beschrijving)
same_descriptions = get_same_descriptions(lib_cluster_info)
count_descriptions = count_descriptions(same_descriptions)


print(same_descriptions)
#print(lib_cluster_info)
#print(np.unique(lib_cluster_info[0]))
