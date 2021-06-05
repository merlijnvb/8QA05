import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re

'''variables'''
verwijderen = ['protein','similar', 'acidic', '4-like', '8-like', '-like', 'ESTs',
                      'like', 'acid-rich', 'adhesion', 'affinity', 'activity', 'alpha-like', 'anion', 'assembly',
                      'association', 'basic', 'candidate', 'carbon','Coenzyme', 'cofactor', 'coiled-coil',
                      'coiled-coil-helix-coiled-coil-helix', 'cold', 'double', 'domain-containing', 'enabled',
                      'fast', 'four', 'glucose', 'half', 'hand', 'inner', 'insert', 'inorganic', 'isoenzyme',
                      'molecule', 'mouse', 'never', 'neighbor', 'nitrogen', 'omega', 'only', 'organic', 'outer', 'paired',
                      'partner', 'region', 'ring', 'slow', 'similarity', 'system', 'very', '3-like', 'beta-like',
                      'coenzyme', 'complex', 'constant', 'component', 'dependent', 'early', 'light', 'long',
                      'protein-like', 'short', '1-like', 'activated', 'group', 'high', 'nine', 'small', 'cell',
                      'chain', 'heavy', 'with', 'acid', 'alpha', 'beta', 'associated', 'containing', 'gamma',
                      'gene', 'inter', 'rich', 'type', 'repeat']
min_frequency = 2
in_nr_clusters = 1
min_length_substring = 4

def data_inlezen(filename): 
    infile = open(filename)
    data = infile.readlines()
    infile.close()
    
    
    lib_data = {}
    if len(data[0].split()) == 1:
        infile_GD = open(filename)             # leest file in
        inlines_GD = infile_GD.read().split()
        infile_GD.close() 
        
        indices = [i for i, x in enumerate(inlines_GD) if x == "\\\\"]  # als item is \\\\, bewaar de index
        lib_data = {}                       # maakt nieuwe dictionary aan
        
        
        line = inlines_GD[: indices[0]]       # elke lijn is het cloneID met zijn beschrijving
        cloneID = int(line[0])                              # geeft cloneID als integer
        
        beschrijving = ""
        for word in range(1, len(line[1:])+1):   
            beschrijving += line[word] + ' '                                    # voegt woord van beschrijving toe aan lijst beschrijving
    
        
        lib_data[cloneID] = beschrijving                              # beschrijving wordt key in dictionary
        
        for i in range(1,len(indices)):   
            line = inlines_GD[indices[i-1]+1: indices[i]]       # elke lijn is het cloneID met zijn beschrijving
            cloneID = int(line[0])                              # geeft cloneID als integer
            
            beschrijving = ""
            for word in range(1, len(line[1:])+1):   
                beschrijving += line[word] + ' '                                    # voegt woord van beschrijving toe aan lijst beschrijving
            
            lib_data[cloneID] = beschrijving                              # beschrijving wordt key in dictionary
    
    else:
        for line in data:
            if line == 'CloneID Familienummer\n':
                pass
            elif len(line.split()[1:]) == 1:                      #als op een regel de ID en cluster nummer staat wordt deze loop gebruikt
                lib_data[int(line.split()[0])] = int(line.split()[1])   #voeg cluster nummer toe aan library, met als key het ID en value een int van clusternummer
            else:                                               #als op een regel de ID en relatieve expressie waarden staan wordt deze gebruikt.
                values = []
                values_line = line.split()
                for value in values_line[1:]:
                    values.append(float(value))
                lib_data[int(values_line[0])] = values  #voeg cluster nummer toe aan library, met als key het ID en value een lijst van floats, met elke float het relatieve expressie waarde
    
    return lib_data

def plot_clusters(expression_data, results):    
    df_data = pd.DataFrame(expression_data)
    days = [1,2,4,7,14,21,45,90] # --> moet geimporteerd worden
    new_cols = {i:days[i] for i in df_data.index}
    df_data = df_data.transpose().rename(columns=new_cols)  
    
    gen_IDs = list(df_data.index.values)
    id_with_cluster = {}
    for ID in gen_IDs:
        id_with_cluster[ID] = results[ID]
    df_data['cluster'] = id_with_cluster.values()
    max_cluster = np.unique(list(id_with_cluster.values())).max()

    df_grouped = df_data.groupby('cluster')
    lower_bound = df_grouped.quantile(q=0.25) # maak eventueel variabel als dit nodig is
    higher_bound = df_grouped.quantile(q=0.75) # maak eventueel variabel als dit nodig is
    middle_bound = df_grouped.quantile(q=0.5) # maak eventueel variabel als dit nodig is
    mean = df_grouped.mean()
    fig, ax=plt.subplots(ncols=1,nrows=max_cluster,figsize=(15,40),sharex=True)
    
    for cluster in mean.index:
        lower_bound.loc[cluster].plot(ax=ax[cluster-1], alpha=1, linewidth=1)
        middle_bound.loc[cluster].plot(ax=ax[cluster-1], alpha=0.5, linewidth=1)
        higher_bound.loc[cluster].plot(ax=ax[cluster-1], alpha=1, linewidth=1)
        mean.loc[cluster].plot(ax=ax[cluster-1], alpha=1, linewidth=3)
        ax[cluster-1].legend(['q=0.25','q=0.5','q=0.75','mean'])
        if cluster == 0:
            cluster = mean.index[-1] + 1
        ax[cluster-1].set_title(f'Cluster: {cluster}')
        ax[cluster-1].set_xticks(days)
        ax[cluster-1].set_xlim(days[0],days[-1])
    
def telwoorden(cluster_data, beschrijvingen_data, verwijderen, min_frequency, in_nr_clusters, min_length_substring):
    cluster_description = {}
    
    for i in range(min(cluster_data.values()),max(cluster_data.values())+1):
        cluster_description[i] = list()                                       # maakt een lijst van de waardes in de range
        
    for ID in cluster_data:
        cluster_description[cluster_data[ID]].append(beschrijvingen_data[ID])             # voegt beschrijvingen toe in de lijst in de dictionary
    
    clusters = list(cluster_description.keys())                                            # maakt lijst van keys in cluster_discription
    descriptions = list(cluster_description.values())                                      # maakt lijst van values in cluster_discription
    lib_substrings = {}                                                     # maakt nieuwe lege dictionary

    '''het volgende stuk (tot de witregel) plaatst alle eerste (of tweede als
    het eerste woord geen woord is, maar één enkele letter/cijfer) woorden als
    key in de dictionary lib_substring'''
    for clust in clusters:
        data_to_check = descriptions[clust-1]                                 
        for desc in data_to_check:
            desc = desc.replace('(', '')
            desc = desc.replace(')', '')
            desc_to_check = re.split(', |_|!| ', desc)                    # splitst de woorden
            
            for bes in desc_to_check:
                if (len(bes) > 3) & (bes not in verwijderen):
                    substring = bes
                    lib_substrings[substring] = ''
            
    for substring in lib_substrings:
        substring_in_clusters = {}   # maakt nieuwe lege dictionary
        for clust in clusters:
            data_to_check_in = descriptions[clust-1]               
            length_substring = len(substring)                    # bepaald lengte substring
            frequency = sum((element[ind:ind+length_substring]).lower() == substring.lower() for element in data_to_check_in for ind,char in enumerate(element)) 
            # telt hoe vaak substrings voorkomen
            if frequency >= min_frequency:
                substring_in_clusters[clust] = frequency
        if len(substring_in_clusters.keys()) == in_nr_clusters:
            lib_substrings[substring] = substring_in_clusters

    keys_to_delete_1 = []                                            # maakt lege lijst om keys te verwijderen    

    for key in lib_substrings:
        if (lib_substrings[key] == "") | (len(key) < min_length_substring):
            keys_to_delete_1.append(key)                            # voegt key toe die verwijdert moet worden aan lijst

    for key in keys_to_delete_1:
        lib_substrings.pop(key)

    return lib_substrings

def plot_familys(data_family, data_expression):  
    
    ID_expr_fam = {}

    for ID in data_family:
        if ID in data_expression.keys():
            expression = data_expression[ID]
            expression.append(int(data_family[ID]))
            ID_expr_fam[ID] = expression

    df_data = pd.DataFrame(ID_expr_fam)
    days = [1,2,4,7,14,21,45,90, 'family'] # --> moet geimporteerd worden
    new_cols = {i:days[i] for i in df_data.index}
    days.pop()
    df_data = df_data.transpose().rename(columns=new_cols).astype(float)
    nr_familys = int(df_data['family'].max())
    
    df_grouped = df_data.groupby('family')
    lower_bound = df_grouped.quantile(q=0.25) # maak eventueel variabel als dit nodig is
    higher_bound = df_grouped.quantile(q=0.75) # maak eventueel variabel als dit nodig is
    middle_bound = df_grouped.quantile(q=0.5) # maak eventueel variabel als dit nodig is
    mean = df_grouped.mean()
    fig, ax=plt.subplots(ncols=1,nrows=len(list(mean.index)),figsize=(30,60),sharex=True)
    plot_nr = 0
    
    for family in mean.index:
        lower_bound.loc[family].plot(ax=ax[plot_nr], alpha=1, linewidth=1)
        middle_bound.loc[family].plot(ax=ax[plot_nr], alpha=0.5, linewidth=1)
        higher_bound.loc[family].plot(ax=ax[plot_nr], alpha=1, linewidth=1)
        mean.loc[family].plot(ax=ax[plot_nr], alpha=1, linewidth=3)
        ax[plot_nr].legend(['q=0.25','q=0.5','q=0.75','mean'])
        if family == 0:
            family = mean.index[-1] + 1
        ax[plot_nr].set_title(f'Family: {family}')
        ax[plot_nr].set_xticks(days)
        ax[plot_nr].set_xlim(days[0],days[-1])
        plot_nr += 1
    
    empty_familys = np.setdiff1d(list(range(1,nr_familys+1)), list(mean.index))

def pie_chart(data_fam, data_clust):
    cluster_fam = {}
    cluster_fam_frequency = {}
    cluster = []
    for ID in data_clust:
        cluster.append(data_clust[ID])
    unique_clusters = np.unique(cluster)
    
    for cluster in unique_clusters:
        cluster_fam_frequency[cluster] = []
    
    for ID in data_fam:
        if ID in data_clust:
            cluster_fam_frequency[data_clust[ID]].append(data_fam[ID])
    
    for cluster in cluster_fam_frequency:
        frequency = {}
        for unique_fam in np.unique(cluster_fam_frequency[cluster]):
            frequency[unique_fam] = cluster_fam_frequency[cluster].count(unique_fam)
        cluster_fam[cluster] = frequency
    
    plot=0
    for cluster in cluster_fam:
        if cluster_fam[cluster] != {}:
            df = pd.DataFrame(cluster_fam[cluster], index=list(range(len(list(cluster_fam[cluster]))))).head(1).transpose()
            df.rename(columns={0: cluster},inplace = True)
            df.plot(kind='pie', subplots=True, title=f'Cluster: {cluster}', autopct='%.0f%%')
            plot += 1
    
    return cluster_fam

lib_cluster_results = data_inlezen('kmca_results.txt')
lib_beschrijvingen = data_inlezen('GenDescription2.txt')
lib_expression = data_inlezen('Voorbeeld_clusterdata.txt')
lib_family = data_inlezen('CloneIdFamily.txt')
plot_clusters(lib_expression, lib_cluster_results)
telwoorden = telwoorden(lib_cluster_results, lib_beschrijvingen, verwijderen, min_frequency, in_nr_clusters, min_length_substring)
plot_familys(lib_family, lib_expression)
pie_chart = pie_chart(lib_family,lib_cluster_results)   
