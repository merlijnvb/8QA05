def TelWoorden(filename):
    
    inFile = open(filename)
    raw_data = inFile.readlines()
    inFile.close()

    data = [[line.strip()] for line in raw_data[1::3]]
    
    desc_freq = {}
    
    for description in data:
        if description not in desc_freq:
            desc_freq[description] = 1
        else:
            desc_freq[description] += 1
            
    desc_freq = {keys: values for keys, values in sorted(desc_freq.items(), key=lambda item: item[1])}
    
    return desc_freq

#%%

def descriptionCloneIDs(gen_description, voorbeeld_cluster_result):
    
    #open gen_description en voorbeeld_cluster_result bestanden
    inFile = open(gen_description)
    raw_data = inFile.readlines()
    inFile.close
    
    voorbeeld_cluster_result = open(voorbeeld_cluster_result)
    raw_cluster_data = voorbeeld_cluster_result.readlines()
    voorbeeld_cluster_result.close()
    
    #maak lijst met de beschrijving en overeenkomende clone IDs
    description = [line.strip()for line in raw_data [1::3]]
    cloneIDs = [line.strip() for line in raw_data in [0::3]]
   
    #maak lege dictionairy genaamd desciptionCloneIDs
    descriptionCloneIDs() = {}
    
    #alle cloneIDs overeenkomend met die beschrijving 
    for item in range(len(descripton)):
        if description[item] not in desciptionCloneIDs:
            descriptionCloneIDs[description[item]] = [cloneIDs[item]]
        else:
            descCloneIDs[description[item]].append[cloneIDs[item]]
   
    #maak lege dictionairy genaamd cluster
    clusterFreq = {}
    
    cluster_data = [items for sublist in [line.strip().splitfor line in raw_cluster_data]for item in sublist]
    
    cloneIDValues = list(descCloneIDs.values())
    cloneIDKeys = list(descCloneIDs.keys())
    
    for description in range(len(cloneIDsValues[description])):
        if cloneIDValues[description][cloneID] in cluster_data:
            cluster_index = cluster_data.index(cloneIDValues[description][cloneID]) + 1
            cluster_value = int(cluster_data[cluster_index])
            clusterFreq[cloneIDKeys[description]][cluster_value] += 1
            
        else: 
            pass
    
    return clusterFreq
        
