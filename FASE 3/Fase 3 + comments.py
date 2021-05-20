# -*- coding: utf-8 -*-
"""
Created on Wed May 19 14:17:20 2021

@author: 20203520
"""


import numpy as np
import re

def lib_res(Filename):
    def openfile(Filename):                     # leest file in
        read_data = open(f'{Filename}')
        data = read_data.read()
        data = data.splitlines()
        read_data.close()
        
        return data
    
    lib = dict()                                # maakt nieuwe dictionary aan
    
    for lines in openfile(Filename):            # splitst de lijnen en voegt de cloneIDs en clusternummers toe aan dictionary
        line = lines.split()
        lib[int(line[0])] = int(line[1])
        
    return lib

def lib_beschrijvingen(Filename):
    infile_GD = open(f'{Filename}')             # leest file in
    inlines_GD = infile_GD.read().split()
    infile_GD.close() 
    
    indices = [i for i, x in enumerate(inlines_GD) if x == "\\\\"]  # als item is \\\\, bewaar de index
    beschrijvingdict = {}                       # maakt nieuwe dictionary aan
    
    
    line = inlines_GD[: indices[0]]       # elke lijn is het cloneID met zijn beschrijving
    cloneID = int(line[0])                              # geeft cloneID als integer
    
    beschrijving = ""
    for word in range(1, len(line[1:])+1):   
        if word == 1:
            beschrijving += line[word]                                    # voegt woord van beschrijving toe aan lijst beschrijving
        
        else:
            beschrijving += " " + line[word]                              # print woorden met spatie
    
    beschrijvingdict[cloneID] = beschrijving                              # beschrijving wordt key in dictionary
    
    for i in range(1,len(indices)):   
        line = inlines_GD[indices[i-1]+1: indices[i]]       # elke lijn is het cloneID met zijn beschrijving
        cloneID = int(line[0])                              # geeft cloneID als integer
        
        beschrijving = ""
        for word in range(1, len(line[1:])+1):   
            if word == 1:
                beschrijving += line[word]                                    # voegt woord van beschrijving toe aan lijst beschrijving

            else:
                beschrijving += " " + line[word]                              # print woorden met spatie
        
        beschrijvingdict[cloneID] = beschrijving                              # beschrijving wordt key in dictionary
        
    return beschrijvingdict

def get_description(input_data, refrence_data):
    lib_refrence = dict()                                                     # maakt nieuwe dictionary aan
    
    for ID in input_data:
        lib_refrence[ID] = refrence_data[ID]                                  # voegt de beschrijvingsdictionary toe aan de nieuwe dictionary
    return lib_refrence

def get_cluster_description(results, refrence_data):
    cluster_discription = dict()                                                # maakt nieuwe dictionary aan
    
    for i in range(min(results.values()),max(results.values())+1):
        cluster_discription[i] = list()                                       # maakt een lijst van de waardes in de range
        
    for ID in results:
        cluster_discription[results[ID]].append(refrence_data[ID])             # voegt beschrijvingen toe in de lijst in de dictionary
    
    return cluster_discription

def multiple_clusters(data, nr_of_clusters, nr_string_in_cluster):
    clusters = list(data.keys())                                            # maakt lijst van keys in cluster_discription
    descriptions = list(data.values())                                      # maakt lijst van values in dluster_discription
    lib_substrings = {}                                                     # maakt nieuwe dictionary
    for clust in clusters:
        data_to_check = descriptions[clust]                                 

        for desc in data_to_check:
            desc_to_check = re.split(', |_|-|!| ', desc)                    # splitst de woorden
            
            # kijkt of nummers/woorden langer zijn dan 1 teken, anders pakt hij het volgende woord/getal
            if (len(desc_to_check[0]) > 1) | len(desc_to_check) == 1:
                substring = desc_to_check[0]
                
            else:
                substring = desc_to_check[1]
                
            lib_substrings[substring] = ''
                
    for substring in lib_substrings:
        substring_in_clusters = {}                                          # maakt nieuwe dictionary
        for clust in clusters:
            #lib_substrings[substring][clust] = sum(substring in list for list in descriptions[clust])
            data_to_check_in = descriptions[clust]               
            length_substring = len(substring)                    # bepaald lengte substring
            count = sum(element[index:index+length_substring] == substring for element in data_to_check_in for index,char in enumerate(element)) 
            
            # telt hoe vaak substrings voorkomen
            if count >= nr_string_in_cluster:
                substring_in_clusters[clust] = count
            #print(substring_in_clusters)
        if len(substring_in_clusters.keys()) == nr_of_clusters:
            lib_substrings[substring] = substring_in_clusters
    
    keys_to_delete = []                                            # maakt lege lijst om keys te verwijderen
    for key in lib_substrings:
        if lib_substrings[key] == "":
            keys_to_delete.append(key)                            # voegt key toe die verwijdert moet worden aan lijst
    
    for key in keys_to_delete:
        lib_substrings.pop(key)                                   # verwijdert key
        
    return lib_substrings


lib_results = lib_res('Voorbeeld_clusterresult.txt')
lib_beschrijving = lib_beschrijvingen('GenDescription2.txt')
lib_bes_toegekend = get_description(lib_results,lib_beschrijving)
lib_cluster_info = get_cluster_description(lib_results, lib_beschrijving)
lib_substrings = multiple_clusters(lib_cluster_info, 3, 2)