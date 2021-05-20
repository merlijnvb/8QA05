import numpy as np
import re

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

def get_description(input_data, refrence_data):
    lib_refrence = dict()
    
    for ID in input_data:
        lib_refrence[ID] = refrence_data[ID]
        
    return lib_refrence

def get_cluster_description(results, refrence_data):
    cluster_discription = dict()
    
    for i in range(min(results.values()),max(results.values())+1):
        cluster_discription[i] = list()
        
    for ID in results:
        cluster_discription[results[ID]].append(refrence_data[ID])
        
    return cluster_discription

def multiple_clusters(data, nr_of_clusters, nr_string_in_cluster):
    clusters = list(data.keys())
    descriptions = list(data.values())
    lib_substrings = {} 
    for clust in clusters:
        data_to_check = descriptions[clust]
        
        for desc in data_to_check:
            desc_to_check = re.split(', |_|-|!| ', desc)
            
            if (len(desc_to_check[0]) > 1) | len(desc_to_check) == 1:
                substring = desc_to_check[0]
            else:
                substring = desc_to_check[1]
                
            lib_substrings[substring] = ''
                
    for substring in lib_substrings:
        substring_in_clusters = {}
        for clust in clusters:
            #lib_substrings[substring][clust] = sum(substring in list for list in descriptions[clust])
            data_to_check_in = descriptions[clust]
            length_substring = len(substring)
            count = sum(element[index:index+length_substring] == substring for element in data_to_check_in for index,char in enumerate(element))
            if count >= nr_string_in_cluster:
                substring_in_clusters[clust] = count
            #print(substring_in_clusters)
        if len(substring_in_clusters.keys()) == nr_of_clusters:
            lib_substrings[substring] = substring_in_clusters
    
    keys_to_delete = []
    for key in lib_substrings:
        if lib_substrings[key] == "":
            keys_to_delete.append(key)
    
    for key in keys_to_delete:
        lib_substrings.pop(key)
        
    return lib_substrings


lib_results = lib_res('Voorbeeld_clusterresult.txt')
lib_beschrijving = lib_beschrijvingen('GenDescription2.txt')
lib_bes_toegekend = get_description(lib_results,lib_beschrijving)
lib_cluster_info = get_cluster_description(lib_results, lib_beschrijving)
lib_substrings = multiple_clusters(lib_cluster_info, 3, 2)