# -*- coding: utf-8 -*-
"""
Created on Mon May 17 16:08:33 2021

@author: 20202139
"""
infile_GD = open('GenDescription2.txt')
gen_beschrijving = infile_GD.readlines()
infile_GD.close()

infile_VC = open('Voorbeeld_clusterresult.txt')
cluster_uitvoer = infile_VC.readlines()
infile_VC.close()

GenDescription = 'Data\GenDescription2.txt'
clusterResultFile = 'Data\Voorbeeld_clusterresult.txt'

def Telwoorden(gen_beschrijving, cluster_uitvoer):
    
    #lijst aanmaken genaams 'items' van alle losse woorden en integers uit 'data'zonder komma's
        
    items = ' '.join([line.strip().lower() for line 
                    in gen_beschrijving[1::3]]).replace('\x01', '')
    
    for ch in [',', '(', '/']:
        if ch in items:
            items = items.replace(ch, ' ')
    items = items.split
    
    cloneIDs = [line.strip() for line in gen_beschrijving[0::3]]
    beschrijvingen = [line.strip().lower() for line in
                      gen_beschrijving [1::3]]

    #lijst aanmaken genaamd 'geen_integers' met uitsluitend woorden uit beschrijvingen
    
    geen_integers = (woord for woord in items if not (woord.isdigit()
                    or woord [0] == '-'and woord[1:].isdigit()))
    
    #lege dictionairy aanmaken genaamd 'woord_freq'
    
    woord_freq ={}
    
    #per woord uit 'geen_integers'de frequentie bepalen en toevoegen aan 'woord_freq'
    
    for woord in geen_integers:
        if woord not in woord_freq:
            woord_freq[woord] = 1
        else:
            woord_freq[woord] += 1
            
    unieke_woorden = list(woord_freq.keys())
    
    beschrijving_cloneIDs = {}
    
    for item in range(len(unieke_woorden)):
        for beschrijving in range(len(beschrijvingen)):
            if '' + unieke_woorden[woord] + ' ' in beschrijvingen[beschrijving]:
                if unieke_woorden[woord] not in beschrijving_cloneIDs:
                    beschrijving_cloneIDs[unieke_woorden[woord]] = [cloneIDs[beschrijving]]
        else:
            beschrijving_cloneIDs[unieke_woorden[woord]].append(cloneIDs[beschrijving])
            
    cloneIDValues = list(beschrijving_cloneIDs.values())
    cloneIDKeys = list(beschrijving_cloneIDs.keys())
     
    cluster_freq = {}
    cluster_data = [item for sublist in [line.strip().split
                    for line in cluster_uitvoer] for item in sublist]
     
    k = (int(max([line for line in cluster_data[1::2]]))+1)
     
    for woord in range(len(cloneIDKeys)):
            
        cluster_freq[cloneIDKeys[woord]] = [0] * k
         
        for cloneID in range(len(cloneIDValues[woord])):
            if cloneIDValues[woord][cloneID] in cluster_data:
                cluster_index = cluster_data.index(cloneIDValues[woord][cloneID]) + 1
                cluster_value = int(cluster_data[cluster_index])
                cluster_freq[cloneIDKeys[woord]][cluster_value] += 1
            else:
                pass
            
    return cluster_freq
 
    
if __name__ == '__main__':
    print(Telwoorden(gen_beschrijving, cluster_uitvoer))