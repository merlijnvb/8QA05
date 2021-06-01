# -*- coding: utf-8 -*-
"""
Created on Tue Jun  1 19:01:20 2021

@author: 20203928
"""

from Bio import Entrez
import numpy as np
import pandas as pd

def lib(FileName):
    def openFile(FileName):
        read_data = open(f'{FileName}')
        data = read_data.read()
        data = data.splitlines()
        read_data.close()
        
        return data
        
    lib = {}
    
    for line in openFile(FileName):
        if line != 'CloneID Familienummer':
            line_list = line.split()
            if len(line_list[1:]) == 1:
                lib[int(line_list[0])] = int(line_list[1])
            else:
                lib[int(line_list[0])] = line_list[1:]
            
    return lib
            
lib_family = lib('..\Data\CloneIdFamily.txt')
lib_clust = lib('..\Data\kmca_results.txt')

cluster = 3

ID_in_cluster = []
for ID in lib_clust :
    if lib_clust[ID] == cluster:
        ID_in_cluster.append(ID)

f = open("../data/accessionnumbers.txt", "r")
lines = f.read().splitlines()

relavant_lines = []
for line in range(len(lines)):
    if line != 0:
        relavant_lines.append(lines[line])

split_line = []
for line in relavant_lines:
    list_line = line.split('\t')
    if '' in list_line:
        list_line.pop()
    split_line.append(list_line)

ID_code = {}
for line in split_line:
    if len(line) == 4:
        ID_code[int(line[0])] = line[-1] 

IDs_in_cluster = []
for ID in lib_clust:
    if lib_clust[ID] == cluster:
       IDs_in_cluster.append(int(ID))

relevant_IDs_code = {}

for ID in IDs_in_cluster:
    if ID in ID_code:
        relevant_IDs_code[ID] = ID_code[ID]

all_locations_gen = []
for ensemblegeneid in list(relevant_IDs_code.values()):
    Entrez.email = 'd.devetzis@student.tue.nl'
    handle = Entrez.efetch(db = 'gene', id=ensemblegeneid) 
    
    lines = handle.readlines()
   
    for index in range(len(lines)):
        if '"Tissue List"' in lines[index]:
            start = index+1
        elif 'label "Category"' in lines[index]:
            stop = index-4
    location_list = lines[start:stop]
    
    
    
    for item in range(len(location_list)):
        if item == 0:
            list_first_item = list(location_list[item])
            if '"' in list_first_item:
                begin_of_locations = list_first_item.index('"')
                locations = location_list[item][begin_of_locations+1:]
                location_list[0] = locations
            for loc in location_list[item].split(';'):
                all_locations_gen.append(loc.strip())
    if all_locations_gen != []:
        if all_locations_gen[-1] == '':
            all_locations_gen.pop()


lower_loc = []
for element in all_locations_gen:
    lower_loc.append(element.lower())
keys = np.unique(lower_loc)
amount_dict = {}

for key in keys:
    amount_dict[key] = lower_loc.count(key)

df = pd.DataFrame(amount_dict, index=list(amount_dict.keys())).head(1).transpose()

df_to_use = df[df['adrenal'] > 5].copy()
#print(df_to_use)

df_to_use.plot(kind='bar', legend=False)

