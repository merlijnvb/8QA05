"""
Course code: 8QA05
group number: 9
members: Dimo Devetzis, Elizabeth Ninh en Femke Schapendonk
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re
from Bio import Entrez
import numpy as np
import pandas as pd
Entrez.email = 'dimodevetzis@gmail.com'

# variables
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




def data_inlezen(filename): 
    infile = open(filename)
    data = infile.readlines()
    infile.close()
    
    
    lib_data = {}
    if len(data[0].split()) == 1:
        infile_GD = open(filename)                                              # read the file
        inlines_GD = infile_GD.read().split()
        infile_GD.close() 
        
        indices = [i for i, x in enumerate(inlines_GD) if x == "\\\\"]          # if item is \\\\, save the index
        lib_data = {}                                                           # creates a new dictionary
        
        
        line = inlines_GD[:indices[0]]                                          # every line is a CloneID with its description
        cloneID = int(line[0])                                                  # makes CloneID an integer
        
        beschrijving = ""
        for word in range(1, len(line[1:])+1):   
            beschrijving += line[word] + ' '                                    # adds a word from the description to the list 'beschrijving'
    
        
        lib_data[cloneID] = beschrijving                                        # description becomes a value in dictionary
        
        for i in range(1,len(indices)):   
            line = inlines_GD[indices[i-1]+1: indices[i]]                       # every line is the CloneID with its description
            cloneID = int(line[0])                                              # gives CloneID as integer
            
            beschrijving = ""
            for word in range(1, len(line[1:])+1):   
                beschrijving += line[word] + ' '                                # adds word from description to list 'beschrijvingen'
            
            lib_data[cloneID] = beschrijving                                    # description becomes value in dictionary
    
    else:
        for line in data:
            if line == 'CloneID Familienummer\n':
                pass
            elif len(line.split()[1:]) == 1:                                    # if  a line consists of the ID and cluster number, this loop is used
                lib_data[int(line.split()[0])] = int(line.split()[1])           # adds cluster number to a dictionary, with as key the ID and value an integer of cluster number
            else:                                                               # if a line has the ID and relative expression, the line is used
                values = []
                values_line = line.split()
                for value in values_line[1:]:
                    values.append(float(value))
                lib_data[int(values_line[0])] = values                          # adds cluster number to a dictionary, with as key the ID and the value a list of floats with the relative expression value
    
    return lib_data

def plot_clusters(expression_data, results):
    pop = []
    for key in expression_data:
        if key not in list(results.keys()):
            pop.append(key)
    for key in pop:
        expression_data.pop(key)
    
    df_data = pd.DataFrame(expression_data)                                     # create a dataframe from the expression data
    days = [1,2,4,7,14,21,45,90]                                                # --> must be imported
    new_cols = {i:days[i] for i in df_data.index}                               # dictionary of the days of the corresponding data points
    df_data = df_data.transpose().rename(columns=new_cols)                      # transpose the dataframe, and put the days as column names
    
    gen_IDs = list(df_data.index.values)                                        # list of all gene IDs in the dataframe   
    id_with_cluster = {}                            
    for ID in gen_IDs:                          
        if ID in results:
            id_with_cluster[ID] = int(results[ID])                              # dictionary with gene ID as key and the cluster number as value 
    # print(id_with_cluster)
    df_data['cluster'] = id_with_cluster.values()                               # adding a new column, which gives for every gene the cluster number 
    max_cluster = len(np.unique(list(id_with_cluster.values())))                # saving the length of the cluster

    df_grouped = df_data.groupby('cluster')                                     # grouping the dataframe by cluster number
    lower_bound = df_grouped.quantile(q=0.25)                                   # creating a new variable if neccessary
    higher_bound = df_grouped.quantile(q=0.7)                                    # creating a new variable if neccessary
    middle_bound = df_grouped.quantile(q=0.5)                                   # creating a new variable if neccessary
    mean = df_grouped.mean()
    fig, ax=plt.subplots(ncols=1,nrows=max_cluster,figsize=(20,20),sharex=True)
    
    for cluster in mean.index:
        lower_bound.loc[cluster].plot(ax=ax[cluster], alpha=1, linewidth=1)
        middle_bound.loc[cluster].plot(ax=ax[cluster], alpha=0.5, linewidth=1)
        higher_bound.loc[cluster].plot(ax=ax[cluster], alpha=1, linewidth=1)
        mean.loc[cluster].plot(ax=ax[cluster], alpha=1, linewidth=3)
        ax[cluster].legend(['q=0.25','q=0.5','q=0.75','mean'])
        ax[cluster].set_title(f'Cluster: {cluster}')
        ax[cluster].set_xticks(days)
        ax[cluster].set_xlim(days[0],days[-1])
    
def telwoorden(cluster_data, beschrijvingen_data, ruwe_data_fase, verwijderen, min_frequency, in_nr_clusters, min_length_substring, alleen_afgevallen = True):
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
        infile = open(ruwe_data_fase)
        data = infile.readlines()
        infile.close()
        
        afgevallen_IDs = []
        for line in data:
            if int(line.split()[0]) not in list(cluster_data.keys()):
                afgevallen_IDs.append(int(line.split()[0]))
        
        afgevallen_desc = []
        for ID in afgevallen_IDs:
            afgevallen_desc.append(beschrijvingen_data[ID])
        
            
            
        cluster_description = {}
        
        for i in range(min(cluster_data.values()),max(cluster_data.values())+1):
            cluster_description[i] = list()                                     # makes a list of the values in the range
           
        for ID in cluster_data:
            cluster_description[cluster_data[ID]].append(beschrijvingen_data[ID]) # adds a description in the list in the dictionary
        clusters = list(cluster_description.keys())                             # makes a list of keys in cluster_discription
        descriptions = list(cluster_description.values())                       # makes a list of values in cluster_discription
        lib_substrings = {}                                                     # makes a new empty dictionary
        '''het volgende stuk (tot de witregel) plaatst alle eerste (of tweede als
        het eerste woord geen woord is, maar één enkele letter/cijfer) woorden als
        key in de dictionary lib_substring'''
        if (in_nr_clusters == 1) & alleen_afgevallen:
            clusters = [-1]
        
        for clust in clusters:
            data_to_check = descriptions[clusters.index(clust)]                                 
            for desc in data_to_check:
                desc = desc.replace('(', '')
                desc = desc.replace(')', '')
                desc_to_check = re.split(', |_|!| ', desc)                      # splits the words on the characters that are located before, after and between | 
                
                for bes in desc_to_check:
                    if (len(bes) > 3) & (bes not in verwijderen):
                        substring = bes
                        lib_substrings[substring] = ''                          # makes a key of the substring
        for substring in lib_substrings:
            substring_in_clusters = {}   # maakt nieuwe lege dictionary
            for clust in clusters:
                data_to_check_in = descriptions[clust-1]               
                length_substring = len(substring)                               # determines the length of the substring
                frequency = sum((element[ind:ind+length_substring]).lower() == substring.lower() for element in data_to_check_in for ind,char in enumerate(element))  # counts the frequency of every substring
                if frequency >= min_frequency:                                  # checks if a substring occurs the minimum amount of times
                    substring_in_clusters[clust] = frequency                    # if the minimum frequency is violated, it is added
            if len(substring_in_clusters.keys()) == in_nr_clusters:             # checks if the substring in the amount of clusters occurs as the user desires
                lib_substrings[substring] = substring_in_clusters
    
        keys_to_delete_1 = []                                                   # makes an empty list to remove keys    
    
        for key in lib_substrings:
            if (lib_substrings[key] == "") | (len(key) < min_length_substring):
                keys_to_delete_1.append(key)                                    # adds key to the list if it must be removed
    
        for key in keys_to_delete_1:
            lib_substrings.pop(key)                                             # the keys that need to be removed are popped
        return lib_substrings
    return 'error'

def check_telwoorden(afgevallen_telwoorden, normal_telwoorden):
    sub_afgevallen = list(afgevallen_telwoorden.keys())
    sub_norm = list(normal_telwoorden.keys())
    for sub in sub_afgevallen:
        if sub in sub_norm:
            afgevallen_telwoorden.pop(sub)
    
    return list(afgevallen_telwoorden.keys())

def plot_familys(data_family, data_expression):  
    
    ID_expr_fam = {}

    for ID in data_family:
        if ID in data_expression.keys():                                        # if the ID is in the expression data, it is used
            expression = data_expression[ID]
            expression.append(int(data_family[ID]))
            ID_expr_fam[ID] = expression                                        # ID as key and the corresponding expression values as value

    df_data = pd.DataFrame(ID_expr_fam)                                         # making a dataframe of the dictionary
    days = [1,2,4,7,14,21,45,90, 'family']                                      # --> needs to be imported
    new_cols = {i:days[i] for i in df_data.index}                               # dictionary for the days corresponding to the data points 
    days.pop()                                                                  # pop 'family' in the list 'days'
    df_data = df_data.transpose().rename(columns=new_cols).astype(float)        # transpose dataframe and rename the column to the days
    nr_familys = int(df_data['family'].max())           #
    
    df_grouped = df_data.groupby('family')
    lower_bound = df_grouped.quantile(q=0.25)                                   # makes a variable if neccessary
    higher_bound = df_grouped.quantile(q=0.75)                                  # makes a variable if neccessary
    middle_bound = df_grouped.quantile(q=0.5)                                   # makes a variable if neccessary
    mean = df_grouped.mean()
    fig, ax=plt.subplots(ncols=1,nrows=len(list(mean.index)),figsize=(30,60),sharex=True)
    plot_nr = 0
    
    for family in mean.index:
        lower_bound.loc[family].plot(ax=ax[plot_nr], alpha=1, linewidth=1)      # plots
        middle_bound.loc[family].plot(ax=ax[plot_nr], alpha=0.5, linewidth=1)   # plots
        higher_bound.loc[family].plot(ax=ax[plot_nr], alpha=1, linewidth=1)     # plots
        mean.loc[family].plot(ax=ax[plot_nr], alpha=1, linewidth=3)             # plots
        ax[plot_nr].legend(['q=0.25','q=0.5','q=0.75','mean'])
        if family == 0:
            family = mean.index[-1] + 1
        ax[plot_nr].set_title(f'Family: {family}')
        ax[plot_nr].set_xticks(days)
        ax[plot_nr].set_xlim(days[0],days[-1])
        plot_nr += 1
    
    empty_familys = np.setdiff1d(list(range(1,nr_familys+1)), list(mean.index)) # save the empty familys as a list

    return empty_familys

def pie_chart(data_fam, data_clust):
    cluster_fam = {}
    cluster_fam_frequency = {}
    cluster = []
    for ID in data_clust:
        cluster.append(data_clust[ID])                                          # append all the clusters to a list
    unique_clusters = np.unique(cluster)                                        # get a list of the clusters, without duplicates
    
    for cluster in unique_clusters:
        cluster_fam_frequency[cluster] = []                                     # set all the clusters as key with empty list as values
    
    for ID in data_fam:
        if ID in data_clust:
            cluster_fam_frequency[data_clust[ID]].append(data_fam[ID])          # append the family number of the ID to the empty list in the right cluster
    
    for cluster in cluster_fam_frequency:
        frequency = {}
        for unique_fam in np.unique(cluster_fam_frequency[cluster]):
            frequency[unique_fam] = cluster_fam_frequency[cluster].count(unique_fam) # counts the number of times the family occurs in the cluster, adds it as a value to a dictionary with as key the family number
        cluster_fam[cluster] = frequency                                        # dictionary with cluster number as key, and a dictionary with the family as key and frequenty as value
    
    plot=0
    for cluster in cluster_fam:
        if cluster_fam[cluster] != {}:
            df = pd.DataFrame(cluster_fam[cluster], index=list(range(len(list(cluster_fam[cluster]))))).head(1).transpose() # create a dataframe with only 1 row
            df.rename(columns={0: cluster},inplace = True)
            df.plot(kind='pie', subplots=True, title=f'Cluster: {cluster}')#, autopct='%.0f%%')       # plot the data as pie-chart
            plot += 1
    
    return cluster_fam 

def file_to_lib(FileName, columns=False):
    '''
        Preconditions:  'path/filename' --> the file you want to open and convert to a dictionary.
        Postconditions: A dictionary:
                        --> Key = first integer 
                        --> Value = floats or integers after the first protein ID integer
                            --> If value after protein ID consist of only 1 element --> return value of dictionary as a integer.
                            --> If value after protein ID constist of more then 1 element --> return value of dictionary as np.array with floats as values of this array.
    '''
    read_data = open(FileName)
    opened_data = read_data.read()
    opened_data = opened_data.splitlines()
    read_data.close()
    lib = {}
    
    for lines in opened_data:
        if 'cloneid' not in lines.lower():
            line = lines.split()
            if len(line[1:]) == 1:
                lib[int(line[0])] = int(line[1])
            else:
                try:
                    lib[int(line[0])] = np.array(line[1:]).astype(float)
                except:
                    lib[int(line[0])] = np.array(line[1:]).astype(str)
        else:
            cols = lines.split()
            
    if columns:
        return lib, cols
    else:
        return lib

lib_family = file_to_lib('CloneIdFamily.txt')
lib_clust = file_to_lib('kmca_results.txt')
lib_codes, columns = file_to_lib('accessionnumbers.txt', True)

def get_location(code='ensemble'):
    df_locations = pd.DataFrame(columns=['location','RPKM','cluster'])
    
    def list_to_string(array, desc=False):
        line = ''
        for element in array:
            line += element.rstrip()
        
        return line
    def list_to_int(array):
        if array != []:
            interger = array[0]
        else:
            interger = 1
        return interger
    def convert(array):
        array = [line[line.find('"',0)+1:line.find('"',-1)-1] for line in array]
        return array 
    
    location = [i for i, element in enumerate(columns[1:]) if code in element][0]
    
    plot = 0
    terms_not_found = []
    for cluster in np.unique(list(lib_clust.values())):
        list_RPKM = []
        terms = [value[location] for index, value in lib_codes.items() if (index in lib_clust.keys()) and (lib_clust[index] == cluster) and (len(value) == 3)]
        terms_not_found.append([index for index, value in lib_codes.items() if (index in lib_clust.keys()) and (lib_clust[index] == cluster) and (len(value) != 3)])
        try:
            handle = Entrez.efetch(db='gene', id=f'{terms} [mus musculus]', retmode = 'html')
            lines = handle.readlines()
            
            
            # functions = [s.rstrip() for i,s in enumerate(lines) if ('enables' in lines[i-1]) and ('anchor' in s)]
            # functions = np.unique(convert(functions), return_counts=True) 
            # processes = [s.rstrip() for i,s in enumerate(lines) if ('involved' in lines[i-1]) and ('anchor' in s)]
            # processes = np.unique(convert(processes), return_counts=True)
            # components = [s.rstrip() for i,s in enumerate(lines) if ('located' in lines[i-1]) and ('anchor' in s)]
            # components = np.unique(convert(components), return_counts=True)
        
            locations = [s.rstrip() for i, s in enumerate(lines) if ('RPKM' in s) and ('expression' in s) and ('}' in lines[i+1])]
            locations += [list_to_string(lines[i:i+[x for x, val in enumerate(lines[i:i+20]) if '}' in val][0]]) for i, s in enumerate(lines) if ('RPKM' in s) and ('expression' in s) and ('}' not in lines[i+1])]
        
            for txt in locations:
                count = 0
                txt = txt.split()
                B = [float(s[:s.find(')')]) for i, s in enumerate(txt) if ('RPKM' in txt[i-1])]
                RPKM = txt[txt.index('expression')+2:]
                #print(RPKM)
                description = [s for i, s in enumerate(RPKM) if ('text' not in s) and ('RPKM' not in s) and ('and' not in s) and ('tissues' not in s) and ('other' not in s) and (s.isdigit() == False)]
                for item in description:
                    if (')' in item) & (count == 0) & (len(description)>2):
                        bound_1 = description.index(item)
                        first_location = ' '.join(description[:bound_1])
                        count += 1
                    elif (')' in item) & (len(description)>2):
                        bound_2 = description.index(item)
                        second_location = ' '.join(description[bound_1+1:bound_2])
                        count += 1
                    
                        
                    
                if count == 1:
                    list_RPKM.append((first_location,B[0]))
                    del first_location
                elif count == 2:
                    list_RPKM.append((first_location,B[0]))
                    list_RPKM.append((second_location,B[1]))
                    del first_location
                    del second_location
            
            sum_lib = {}
            
            for loc in list_RPKM:
                keys = sum_lib.keys()
                for key in keys:
                    if loc[0].lower() in key:
                        sum_lib[key] = sum_lib[key]+loc[1]
                else:
                     sum_lib[loc[0].lower()] = loc[1]
            #print(sum_lib)
            
        except:
            sum_lib = {}
            for time in range(2):
                if time == 0:
                    terms_to_use = terms[:int(0.5*len(terms))]
                else:
                    terms_to_use = terms[int(0.5*len(terms)):]
                handle = Entrez.efetch(db='gene', id=f'{terms_to_use} [mus musculus]', retmode = 'html')
                lines = handle.readlines()
            
                
                # functions = [s.rstrip() for i,s in enumerate(lines) if ('enables' in lines[i-1]) and ('anchor' in s)]
                # functions = np.unique(convert(functions), return_counts=True) 
                # processes = [s.rstrip() for i,s in enumerate(lines) if ('involved' in lines[i-1]) and ('anchor' in s)]
                # processes = np.unique(convert(processes), return_counts=True)
                # components = [s.rstrip() for i,s in enumerate(lines) if ('located' in lines[i-1]) and ('anchor' in s)]
                # components = np.unique(convert(components), return_counts=True)
            
                locations = [s.rstrip() for i, s in enumerate(lines) if ('RPKM' in s) and ('expression' in s) and ('}' in lines[i+1])]
                locations += [list_to_string(lines[i:i+[x for x, val in enumerate(lines[i:i+20]) if '}' in val][0]]) for i, s in enumerate(lines) if ('RPKM' in s) and ('expression' in s) and ('}' not in lines[i+1])]
            
                for txt in locations:
                    count = 0
                    txt = txt.split()
                    B = [float(s[:s.find(')')]) for i, s in enumerate(txt) if ('RPKM' in txt[i-1])]
                    RPKM = txt[txt.index('expression')+2:]
                    description = [s for i, s in enumerate(RPKM) if ('text' not in s) and ('RPKM' not in s) and ('and' not in s) and ('tissues' not in s) and ('other' not in s) and (s.isdigit() == False)]
                    for item in description:
                        if (')' in item) & (count == 0) & (len(description)>2):
                            bound_1 = description.index(item)
                            first_location = ' '.join(description[:bound_1])
                            count += 1
                        elif (')' in item) & (len(description)>2):
                            bound_2 = description.index(item)
                            second_location = ' '.join(description[bound_1+1:bound_2])
                            count += 1
                        
                            
                        
                    if count == 1:
                        list_RPKM.append((first_location,B[0]))
                        del first_location
                    elif count == 2:
                        list_RPKM.append((first_location,B[0]))
                        list_RPKM.append((second_location,B[1]))
                        del first_location
                        del second_location
                
                for loc in list_RPKM:
                    keys = sum_lib.keys()
                    for key in keys:
                        if loc[0].lower() in key:
                            sum_lib[key] = sum_lib[key]+loc[1]
                    else:
                         sum_lib[loc[0].lower()] = loc[1]



        df = pd.DataFrame(sum_lib, index=list(sum_lib.keys())).head(1).transpose()
        df_to_use = df[df[list(df)[0]]>150].copy()            
        for i in range(len(df_to_use)):                    
            df_locations = df_locations.append({'location': df_to_use.index[i],'RPKM': df_to_use.values[i][0] ,'cluster': cluster}, ignore_index=True)

        df_to_use.plot(kind='bar', legend=False, title=cluster)
        plot += 1

        
    print(df_locations.set_index(['cluster','location']))
    return terms_not_found
def main(f_clus_res, f_desc, f_exprs, f_fam, min_frequency, in_nr_clusters, min_length_substring, fase1_data):
    lib_cluster_results = data_inlezen(f_clus_res)
    lib_beschrijvingen = data_inlezen(f_desc)
    lib_expression = data_inlezen(f_exprs)
    
    if input('Would you like to know the gene distribution of the clusters? [y]/n \n') == 'y':
        plot_clusters(lib_expression, lib_cluster_results)
   
    question = input(f'How many clusters would you like to check for similar descriptions? \nChoose a number > 0 & <= {len(np.unique(list(lib_cluster_results.values())))} \nIf you don\'t want to check, type \"n\" \n')
    if question != "n":
        in_nr_clusters = question
        in_nr_clusters = int(in_nr_clusters)   
        telwoorden_res = telwoorden(lib_cluster_results, lib_beschrijvingen, f_exprs, verwijderen, min_frequency, in_nr_clusters, min_length_substring, alleen_afgevallen=False)
        if telwoorden_res != 'error':
            print(telwoorden_res)
        else:
            return
    
    if input('Would you like to know which descriptions appear almost exclusively in the unclustered genes? [y]/n? \n') == 'y':
        
        afgevallen = telwoorden(lib_cluster_results, lib_beschrijvingen, f_exprs, verwijderen, 1, 1, min_length_substring, alleen_afgevallen=True)
        normal = telwoorden(lib_cluster_results, lib_beschrijvingen, f_exprs, verwijderen, 1, 1, min_length_substring, alleen_afgevallen=False)
        print(f'The descriptions which appear almost exclusively in the unclustered genes are: {check_telwoorden(afgevallen, normal)}')
   
    if input('Would you like to know the changes over time of the families? [y]/n? \n') == 'y':
        lib_family = data_inlezen(f_fam)
        empty_familys = plot_familys(lib_family, lib_expression)
        print(f'The following families have not been assigned to a cluster: {empty_familys}')
        
        if input('Would you like the distribution over the clusters in pie-charts? [y]/n \n') =='y':
            empty_familys = pie_chart(lib_family, lib_cluster_results)

            
    elif input('Would you like the distribution in pie-charts? [y]/n \n')=='y':
        empty_familys = pie_chart(lib_family, lib_cluster_results)
    
    
    if input('Would you like to know the location where the genes of the clusters come to expression? [y]/n\n') == 'y':
        terms_without_ensemblegeneid = get_location()
        if input("Do you want to know the gen ID's for wich no ensemblegeneid is available? [y]/n \n") == 'y':
            for cluster in range(len(terms_without_ensemblegeneid)):
                print(f"In cluster {cluster}, the following ensemblegeneid's where not available: \n\n{terms_without_ensemblegeneid[cluster]} \n\n\n")

min_frequency = 2
in_nr_clusters = 1
min_length_substring = 4
f_clus_res = 'kmca_results.txt'
f_desc = 'GenDescription2.txt'
f_exprs = 'filterd.txt'
f_fam = 'CloneIdFamily.txt'

main(f_clus_res, f_desc, f_exprs, f_fam, min_frequency, in_nr_clusters, min_length_substring, f_exprs)
