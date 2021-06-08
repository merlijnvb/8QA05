"""
Course code: 8QA05
group number: 9
members: Dimo Devetzis, Elizabeth Ninh en Femke Schapendonk
"""

# libraries
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re


def data_inlezen(filename): 
    infile = open(filename)
    data = infile.readlines()
    infile.close()
    
    
    lib_data = {}
    if len(data[0].split()) == 1:
        infile_GD = open(filename)                                              # leest file in
        inlines_GD = infile_GD.read().split()
        infile_GD.close() 
        
        indices = [i for i, x in enumerate(inlines_GD) if x == "\\\\"]          # als item is \\\\, bewaar de index
        lib_data = {}                                                           # maakt nieuwe dictionary aan
        
        
        line = inlines_GD[:indices[0]]                                          # elke lijn is het cloneID met zijn beschrijving
        cloneID = int(line[0])                                                  # geeft cloneID als integer
        
        beschrijving = ""
        for word in range(1, len(line[1:])+1):   
            beschrijving += line[word] + ' '                                    # voegt woord van beschrijving toe aan lijst beschrijving
    
        
        lib_data[cloneID] = beschrijving                                        # beschrijving wordt key in dictionary
        
        for i in range(1,len(indices)):   
            line = inlines_GD[indices[i-1]+1: indices[i]]                       # elke lijn is het cloneID met zijn beschrijving
            cloneID = int(line[0])                                              # geeft cloneID als integer
            
            beschrijving = ""
            for word in range(1, len(line[1:])+1):   
                beschrijving += line[word] + ' '                                # voegt woord van beschrijving toe aan lijst beschrijving
            
            lib_data[cloneID] = beschrijving                                    # beschrijving wordt key in dictionary
    
    else:
        for line in data:
            if line == 'CloneID Familienummer\n':
                pass
            elif len(line.split()[1:]) == 1:                                    # als op een regel de ID en cluster nummer staat wordt deze loop gebruikt
                lib_data[int(line.split()[0])] = int(line.split()[1])           # voeg cluster nummer toe aan library, met als key het ID en value een int van clusternummer
            else:                                                               # als op een regel de ID en relatieve expressie waarden staan wordt deze gebruikt.
                values = []
                values_line = line.split()
                for value in values_line[1:]:
                    values.append(float(value))
                lib_data[int(values_line[0])] = values                          # voeg cluster nummer toe aan library, met als key het ID en value een lijst van floats, met elke float het relatieve expressie waarde
    
    return lib_data

def plot_clusters(expression_data, results):    
    df_data = pd.DataFrame(expression_data)                                     # create a dataframe from the expression data
    days = [1,2,4,7,14,21,45,90]                                                # --> moet geimporteerd worden
    new_cols = {i:days[i] for i in df_data.index}                               # dictionary voor de dagen die horen bij de data punten
    df_data = df_data.transpose().rename(columns=new_cols)                      # transpose the dataframe, and put the days as column names
    
    gen_IDs = list(df_data.index.values)                                        # lijst van alle gen ID's in het dataframe   
    id_with_cluster = {}                            
    for ID in gen_IDs:                          
        id_with_cluster[ID] = int(results[ID])                                  # library met gen ID als key en het cluster nummer als value 
    df_data['cluster'] = id_with_cluster.values()                               # een nieuwe column toevoegen, die voor elk gen het cluster nummer weer geeft. 
    max_cluster = np.unique(list(id_with_cluster.values())).max()               # het grootste cluster nummer opslaan

    df_grouped = df_data.groupby('cluster')                                     # groepeer dataframe bij cluster nummer
    lower_bound = df_grouped.quantile(q=0.25)                                   # maak eventueel variabel als dit nodig is
    higher_bound = df_grouped.quantile(q=0.75)                                  # maak eventueel variabel als dit nodig is
    middle_bound = df_grouped.quantile(q=0.5)                                   # maak eventueel variabel als dit nodig is
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
    
def telwoorden(cluster_data, beschrijvingen_data, ruwe_data_fase, verwijderen, min_frequency, in_nr_clusters, min_length_substring, alleen_afgevallen = True):    
    cluster_description = {}
    if (in_nr_clusters <= 0) and (min_frequency > 0):
        print(f'ERROR: in_nr_clusters moet een getal zijn groter dan 0, is nu: {in_nr_clusters}')
        pass
    elif (in_nr_clusters >= 0) and (min_frequency <= 0):
        print(f'ERROR: min_frequency moet een getal zijn groter dan 0, is nu: {min_frequency}')
        pass
    elif (in_nr_clusters <= 0) and (min_frequency <= 0):
        print(f'ERROR: min_frequency en in_nr_clusters moeten groter zijn  dan 0, \nmin_frequency is nu: {min_frequency} \nin_nr_cluster is nu: {in_nr_clusters}')
        pass
    else:
        for i in range(min(cluster_data.values()),max(cluster_data.values())+1):
            cluster_description[i] = list()                                         # maakt een lijst van de waardes in de range
           
        for ID in cluster_data:
            cluster_description[cluster_data[ID]].append(beschrijvingen_data[ID])   # voegt beschrijvingen toe in de lijst in de dictionary
        clusters = list(cluster_description.keys())                                 # maakt lijst van keys in cluster_discription
        descriptions = list(cluster_description.values())                           # maakt lijst van values in cluster_discription
        
        if (in_nr_clusters == 1) & alleen_afgevallen:                               # als er gevraagd wordt naar de beschrijving van de afgevallen genen, wordt het bestand met de ongefilterde data van fase 1 gebruikt
            infile = open(ruwe_data_fase)
            data = infile.readlines()
            infile.close()
            
            afgevallen_IDs = []
            for line in data:
                if int(line.split()[0]) not in list(cluster_data.keys()):           # als het ID niet voorkomt in de geclusterde ID's, wordt het ID toegevoegd aan de lijst afgevallen_IDs
                    afgevallen_IDs.append(int(line.split()[0]))
            
            afgevallen_desc = []
            for ID in afgevallen_IDs:
                afgevallen_desc.append(beschrijvingen_data[ID])                     # lijst met beschrijvingen van de afgevallen genen
                
            clusters = [-1]                                                         # lijst met afgevallen genen, voor het bekijken van hun beschrijvingen
            
        '''het volgende stuk (tot de witregel) plaatst alle eerste (of tweede als
        het eerste woord geen woord is, maar één enkele letter/cijfer) woorden als
        key in de dictionary lib_substring'''
        
        lib_substrings = {}                                                         # maakt nieuwe lege dictionary
    
        for clust in clusters:
            data_to_check = descriptions[clusters.index(clust)]                                 
            for desc in data_to_check:
                desc = desc.replace('(', '')
                desc = desc.replace(')', '')
                desc_to_check = re.split(', |_|!| ', desc)                          # splitst de woorden op de caracters die voor, na en tussen | staan
                
                for bes in desc_to_check:
                    if (len(bes) > 3) & (bes not in verwijderen):
                        substring = bes
                        lib_substrings[substring] = ''                              # maakt een key van de substrings
        for substring in lib_substrings:
            substring_in_clusters = {}                                              # maakt nieuwe lege dictionary
            for clust in clusters:
                data_to_check_in = descriptions[clust-1]               
                length_substring = len(substring)                                   # bepaald lengte substring
                frequency = sum((element[ind:ind+length_substring]).lower() == substring.lower() for element in data_to_check_in for ind,char in enumerate(element))  # telt hoe vaak substrings voorkomen
                if frequency >= min_frequency:                                      # kijkt of een substring minimaal het aantal keer voorkomt als dat we willen
                    substring_in_clusters[clust] = frequency                        # als de minimale frequency is overschreden, wordt die toegevoegd
            if len(substring_in_clusters.keys()) == in_nr_clusters:                 # checkt of de substring in het aantal clusters voorkomt dat je wilt
                lib_substrings[substring] = substring_in_clusters
    
        keys_to_delete_1 = []                                                       # maakt lege lijst om keys te verwijderen    
    
        for key in lib_substrings:
            if (lib_substrings[key] == "") | (len(key) < min_length_substring):
                keys_to_delete_1.append(key)                                        # voegt key toe die verwijdert moet worden aan lijst
    
        for key in keys_to_delete_1:
            lib_substrings.pop(key)                                                 # de keys die verwijderd moeten worden worden gepopt
        return lib_substrings
    return 'error'

def check_telwoorden(afgevallen_telwoorden, normal_telwoorden):
    ''' Vergelijkt de beschrijvingen van de wel-geclusterde genen met die van de niet-geclusterde genen, 
    en returnt een lijst met woorden die vrijwel alleen voorkomen in de niet-geclusterde genen'''
    sub_afgevallen = list(afgevallen_telwoorden.keys())
    sub_norm = list(normal_telwoorden.keys())
    for sub in sub_afgevallen:
        if sub in sub_norm:                                                     # als de beschrijving wel voorkomt in de wel-geclusterde genen, haal deze beschrijving dan uit de lijst afgevallen_telwoorden
            afgevallen_telwoorden.pop(sub)
    
    return list(afgevallen_telwoorden.keys())

def plot_familys(data_family, data_expression):  
    
    ID_expr_fam = {}

    for ID in data_family:
        if ID in data_expression.keys():                                        # als het ID in de expressie data zit dan wordt die gebruikt
            expression = data_expression[ID]
            expression.append(int(data_family[ID]))
            ID_expr_fam[ID] = expression                                        # ID als key en de expressiewaardes voor dat ID als value

    df_data = pd.DataFrame(ID_expr_fam)                                         # dataframe maken van de library
    days = [1,2,4,7,14,21,45,90, 'family']                                      # --> moet geimporteerd worden
    new_cols = {i:days[i] for i in df_data.index}                               # dictionary voor de dagen die horen bij de data punten 
    days.pop()                                                                  # pop 'family' van de lijst days
    df_data = df_data.transpose().rename(columns=new_cols).astype(float)        # transpose dataframe en hernoem de colum naar de dagen
    nr_familys = int(df_data['family'].max())           #
    
    df_grouped = df_data.groupby('family')
    lower_bound = df_grouped.quantile(q=0.25)                                   # maak eventueel variabel als dit nodig is
    higher_bound = df_grouped.quantile(q=0.75)                                  # maak eventueel variabel als dit nodig is
    middle_bound = df_grouped.quantile(q=0.5)                                   # maak eventueel variabel als dit nodig is
    mean = df_grouped.mean()
    fig, ax=plt.subplots(ncols=1,nrows=len(list(mean.index)),figsize=(30,60),sharex=True)
    plot_nr = 0
    
    for family in mean.index:
        lower_bound.loc[family].plot(ax=ax[plot_nr], alpha=1, linewidth=1)      # plotten
        middle_bound.loc[family].plot(ax=ax[plot_nr], alpha=0.5, linewidth=1)   # plotten
        higher_bound.loc[family].plot(ax=ax[plot_nr], alpha=1, linewidth=1)     # plotten
        mean.loc[family].plot(ax=ax[plot_nr], alpha=1, linewidth=3)             # plotten
        ax[plot_nr].legend(['q=0.25','q=0.5','q=0.75','mean'])
        if family == 0:
            family = mean.index[-1] + 1
        ax[plot_nr].set_title(f'Family: {family}')
        ax[plot_nr].set_xticks(days)
        ax[plot_nr].set_xlim(days[0],days[-1])
        plot_nr += 1
    
    empty_familys = np.setdiff1d(list(range(1,nr_familys+1)), list(mean.index)) # sla de lege families op als een lijst 

    return empty_familys

def pie_chart(data_fam, data_clust):
    cluster_fam = {}
    cluster_fam_frequency = {}
    cluster = []
    for ID in data_clust:
        cluster.append(data_clust[ID])                                          # voeg alle clusters toe aan een lijst 
    unique_clusters = np.unique(cluster)                                        # maak een lijst van clusters zonder duplicates
    
    for cluster in unique_clusters:
        cluster_fam_frequency[cluster] = []                                     # maak een dictionary met als key het cluster en een lege lijst als value 
    
    for ID in data_fam:
        if ID in data_clust:
            cluster_fam_frequency[data_clust[ID]].append(data_fam[ID])          # voeg het familienummer van het ID toe aan de lege lijst in het juiste cluster
    
    for cluster in cluster_fam_frequency:
        frequency = {}
        for unique_fam in np.unique(cluster_fam_frequency[cluster]):
            frequency[unique_fam] = cluster_fam_frequency[cluster].count(unique_fam) # telt hoevaak family in dat cluster voorkomt, voegd dit toe (als value) aan lib, met als key het familie nummer
        cluster_fam[cluster] = frequency                                        # dictionary met cluster nummer als key, en een dictionary, met de familie als key en frequentie als value, als value
    
    plot = 0
    for cluster in cluster_fam:
        if cluster_fam[cluster] != {}:
            df = pd.DataFrame(cluster_fam[cluster], index=list(range(len(list(cluster_fam[cluster]))))).head(1).transpose() # maak een dataframe met maar één rij
            df.rename(columns={0: cluster},inplace = True)
            df.plot(kind='pie', subplots=True, title=f'Cluster: {cluster}', autopct='%.0f%%') # plot de data als pie-chart
            plot += 1
    
    return cluster_fam 

'''vrije interpretatie moet hier nog komen'''

def main(f_clus_res, f_desc, f_exprs, f_fam, min_frequency, in_nr_clusters, min_length_substring, verwijderen,fase1_data):
    lib_cluster_results = data_inlezen(f_clus_res)
    lib_beschrijvingen = data_inlezen(f_desc)
    lib_expression = data_inlezen(f_exprs)
    
    if input('Would you like to know which descriptions appear almost exlcusively in the unclustered genes? [y]/n? \n') == 'y': 
        afgevallen = telwoorden(lib_cluster_results, lib_beschrijvingen, fase1_data, verwijderen, min_frequency, in_nr_clusters, min_length_substring, alleen_afgevallen=True)
        normal = telwoorden(lib_cluster_results, lib_beschrijvingen, fase1_data, verwijderen, min_frequency, in_nr_clusters, min_length_substring, alleen_afgevallen=False)
        print(f'de beschrijvingen die vrijwel alleen voorkomen in de niet geclusterde genen zijn: {check_telwoorden(afgevallen, normal)}')
        
    in_nr_clusters = input(f'How many clusters would you like to check for similar discriptions? \nChoose a number > 0 & <= {max(list(lib_cluster_results.values()))} \nIf you don\'t want to check, type \"n\" \n')
    if in_nr_clusters != "n":
        in_nr_clusters = int(in_nr_clusters)   
        telwoorden_res = telwoorden(lib_cluster_results, lib_beschrijvingen, fase1_data, verwijderen, min_frequency, in_nr_clusters, min_length_substring, alleen_afgevallen=False)
        if telwoorden_res != 'error':
            print(telwoorden_res)
        else:
            return
        
    if input("Would you like the plots of phase 3? [y]/n \n")=="y":
        if input('Would you like to know the gene distribution of the clusters? [y]/n \n') == 'y':
            plot_clusters(lib_expression, lib_cluster_results)
       
        if input('Would you like to know the family distribution of the clusters? [y]/n \n') == 'y':
            lib_family = data_inlezen(f_fam)
    
            if input('Would you like them in a lineplot? [y]/n \n')=='y':
                empty_familys = plot_familys(lib_family, lib_expression)
                print(f'The following families have not been assigned to a cluster: {empty_familys}')
                
                if input('Would you like to see the family distribution over the cluster in pie charts? [y]/n \n') =='y':
                    empty_familys = pie_chart(lib_family, lib_cluster_results)
    
                
            elif input('Would you like the distribution in pie-charts? [y]/n \n')=='y':
                empty_familys = pie_chart(lib_family, lib_cluster_results)

    return 
