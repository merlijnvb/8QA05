# -*- coding: utf-8 -*-
"""
Created on Mon May 17 15:17:02 2021

@author: 20202876
"""

GenDescription = 'Data\GenDescription.txt'
clusterResultFile = 'Data\Voorbeeld_clusterresult.txt'

def Telwoorden(gen_beschrijving, cluster_uitvoer):
    
    #inlezen van 'gen_beschrijving'
    #rauwe data uit 'gen_beschrijving' toewijzen aan variabele 'gen beschrijving'
    #rauwe data uit 'cluster_uitvoer' toewijzen aan variabele 'cluster_uitvoer'
    
    with open(gen_beschrijving) as gen_beschrijving:
        gen_beschrijving = gen_beschrijving.readlines()
        
    with open(cluster_uitvoer) as cluster_uitvoer
        cluster_uitvoer = cluster_uitvoer.readlines()
        
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
    
    geen_integers = [ woord for woord in items if not (wood.isdigit()
                    or woord [0] == '-'and woord[1:].isdigit())]
    
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
    
     for item in range(len(unieke woorden)):
        for beschrijving in range(len(beschrijvingen)):
        if '' + unieke_woorden[woord] + ' ' in beschrijvingen[beschrijving]:
            if unieke_woorden[woord] not in beschrijving_cloneIDs:
            beschrijving_cloneIDs[unieke_woorden[woord]] = [cloneIDs[beschrijving]]
        else:
            beschrijving_cloneIDs[unieke_woorden[woord]].append(cloneIDs[beschrijving])
            
     cloneIDValues = list(beschrijving_cloneID.values())
     cloneIDKeys = list(beschrijving_cloneID.keys())
     
     cluster_freq = {}
     cluster_data = [item for sublist in [line.strip().split
                    for line in cluster_uitvoer] for item in sublist]
     
     k = (int(max([line for line in cluster_data[1::2]]))+1)
     
     for woord in range(len(cloneIDs_keys)):
            
         cluste_freq[cloneID_keys[woord]] = [0] * k
         
         for cloneID in range(len(cloneID_values[woord])):
             if cloneID_values[woord][cloneID] in cluster_data:
                 cluster_index = cluster_data.index(cloneID_values[woord][cloneID]) + 1
                 cluser_value = int(cluster_data[cluster_index])
                 cluster_freq[cloneID_keys[woord]][cluster_value] += 1
            else:
                pass
            
     return cluster_freq
 
    
if _name_ == '_main_':
    print(TelWoorden(GenDescription, clusterResultFile))
     
                 
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
                      


