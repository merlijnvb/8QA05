import numpy as np
import re

def lib_res(Filename):
    '''
    preconditions:
        - Filename = name of the file with the cluster results
    postconditions: 
        - returns a library with the genID as key and cluster as value
    '''
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
    '''
    preconditions:
        - Filename = name of the file with descriptions
    postconditions:
        library with genID as key and description of gen as value'''
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

#def get_description(input_data, refrence_data):
    # lib_refrence = dict()                                                     # maakt nieuwe dictionary aan
    
    # for ID in input_data:
    #     lib_refrence[ID] = refrence_data[ID]                                  # voegt de beschrijvingsdictionary toe aan de nieuwe dictionary
    # return lib_refrence

def get_cluster_description(results, refrence_data):
    cluster_discription = dict()                                                # maakt nieuwe dictionary aan
    
    for i in range(min(results.values()),max(results.values())+1):
        cluster_discription[i] = list()                                       # maakt een lijst van de waardes in de range
        
    for ID in results:
        cluster_discription[results[ID]].append(refrence_data[ID])             # voegt beschrijvingen toe in de lijst in de dictionary
    
    return cluster_discription

def telwoorden(data, nr_of_clusters, nr_string_in_cluster):
    '''
    preconditions: 
        - data --> library:
                    key   = cluster
                    value = beschrijving
        - nr_of_clusters --> in hoeveel clusters de substring voor moet komen
        - nr_string_in_cluster --> hoevaak de substring minimaal in het 
                                    cluster voor moet komen
                                    
    postconditions:
        return een library (lib_substrings):
            key   --> substring
            value --> library:
                key   --> cluster
                value --> frequency
    taak:
        Deze functie bepaald hoevaak een substring voorkomt in de clusters
    '''
    clusters = list(data.keys())                                            # maakt lijst van keys in cluster_discription
    descriptions = list(data.values())                                      # maakt lijst van values in cluster_discription
    lib_substrings = {}                                                     # maakt nieuwe lege dictionary
    
    '''het volgende stuk (tot de witregel) plaatst alle eerste (of tweede als
    het eerste woord geen woord is, maar één enkele letter/cijfer) woorden als
    key in de dictionary lib_substring'''
    for clust in clusters:
        data_to_check = descriptions[clust]                                 
        for desc in data_to_check:
            # print(desc)
            desc = desc.replace('(', '')
            desc = desc.replace(')', '')
            desc_to_check = re.split(', |_|!| ', desc)                    # splitst de woorden
            
            for bes in desc_to_check:
                if len(bes) > 3:
                    substring = bes
                lib_substrings[substring] = ''
            
    for substring in lib_substrings:
        substring_in_clusters = {}   # maakt nieuwe lege dictionary
        for clust in clusters:
            data_to_check_in = descriptions[clust]               
            length_substring = len(substring)                    # bepaald lengte substring
            frequency = sum((element[ind:ind+length_substring]).lower() == substring.lower() for element in data_to_check_in for ind,char in enumerate(element)) 
            # telt hoe vaak substrings voorkomen
            if frequency >= nr_string_in_cluster:
                substring_in_clusters[clust] = frequency
        if len(substring_in_clusters.keys()) == nr_of_clusters:
            lib_substrings[substring] = substring_in_clusters
    keys_to_delete = []                                            # maakt lege lijst om keys te verwijderen
    for key in lib_substrings:
        if (lib_substrings[key] == "") | (len(key) < 3):
            keys_to_delete.append(key)                            # voegt key toe die verwijdert moet worden aan lijst
    for key in keys_to_delete:
        lib_substrings.pop(key)                                   # verwijdert key
        
    return lib_substrings

def fucntie_uitvoeren(Filename_res, Filename_bes, clusters, in_cluster):
    '''
    preconditions: 
        - Filename_res = name of the file with the cluster results
        - Filename_bes = name of the file with the descriptions of each gen
        - clusters = substring in how many of clusters
        - in_cluster = minimum frequency of a substring in a cluster
    postconditions:
        - returns the result of the function multiple_clusters 
    taak:
        - alle functies aanroepen
    '''   
    lib_results = lib_res(Filename_res)
    lib_beschrijving = lib_beschrijvingen(Filename_bes)
    lib_cluster_info = get_cluster_description(lib_results, lib_beschrijving)
    lib_substrings = telwoorden(lib_cluster_info, clusters, in_cluster)
    
    return lib_substrings, lib_beschrijving

lib_substrings, lib_beschrijving = fucntie_uitvoeren('Voorbeeld_clusterresult.txt', 'GenDescription2.txt', 6,1)
print(lib_substrings)
