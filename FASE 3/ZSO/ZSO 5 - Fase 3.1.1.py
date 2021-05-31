# -*- coding: utf-8 -*-
"""
Created on Tue May 18 11:57:19 2021

@author: 20202139
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
    for word in line[1:]:                   
        beschrijving += " " + word                              # print woorden met spatie
    
    beschrijvingdict[cloneID] = beschrijving
    
    for i in range(1,len(indices)):   
        line = inlines_GD[indices[i-1]+1: indices[i]]       # elke lijn is het cloneID met zijn beschrijving
        cloneID = int(line[0])                              # geeft cloneID als integer
        
        beschrijving = ""
        for word in line[1:]:                   
            beschrijving += " " + word                              # print woorden met spatie
        
        beschrijvingdict[cloneID] = beschrijving
        
    return beschrijvingdict
    
lib_results = lib_res('Voorbeeld_clusterresult.txt')
lib_beschrijving = lib_beschrijvingen('GenDescription2.txt')

def get_descrision(input_data, refrence_data):
    lib_refrence = dict()
    
    for ID in input_data:
        lib_refrence[ID] = refrence_data[ID]
        
    return lib_refrence

lib_bes_toegekend = get_descrision(lib_results,lib_beschrijving)

def get_cluster_descrision(results, refrence_data):
    cluster_discription = dict()
    
    for i in range(min(results.values()),max(results.values())+1):
        cluster_discription[i] = list()
        
    for ID in results:
        cluster_discription[results[ID]].append(refrence_data[ID])
        
    return cluster_discription
    
    
lib_cluster_info = get_cluster_descrision(lib_results, lib_beschrijving)
print(lib_cluster_info)
#print(np.unique(lib_cluster_info[0]))
