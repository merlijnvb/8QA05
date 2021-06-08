from Bio import Entrez
import numpy as np
import time
Entrez.email = 'dimodevetzis@gmail.com' 

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

lib_family = file_to_lib('Data\CloneIdFamily.txt')
lib_clust = file_to_lib('Data\kmca_results.txt')
lib_codes, columns = file_to_lib('Data/accessionnumbers.txt', True)

def get_location(code='ensemble'):    
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
    
    for cluster in [4]:#'''np.unique(list(lib_clust.values())):'''
        terms = [value[location] for index, value in lib_codes.items() if (index in lib_clust.keys()) and (lib_clust[index] == cluster) and (len(value) == 3)]
        
        terms_not_found = [index for index, value in lib_codes.items() if (index in lib_clust.keys()) and (lib_clust[index] == cluster) and (len(value) != 3)]
        #print(cluster, terms_not_found)
        try:
            handle = Entrez.efetch(db='gene', id=f'{terms} [mus musculus]', retmode = 'html')
            lines = handle.readlines()
            
            
            functions = [s.rstrip() for i,s in enumerate(lines) if ('enables' in lines[i-1]) and ('anchor' in s)]
            functions = np.unique(convert(functions), return_counts=True) 
            processes = [s.rstrip() for i,s in enumerate(lines) if ('involved' in lines[i-1]) and ('anchor' in s)]
            processes = np.unique(convert(processes), return_counts=True)
            components = [s.rstrip() for i,s in enumerate(lines) if ('located' in lines[i-1]) and ('anchor' in s)]
            components = np.unique(convert(components), return_counts=True)
        
            locations = [s.rstrip() for i, s in enumerate(lines) if ('RPKM' in s) and ('expression' in s) and ('}' in lines[i+1])]
            locations += [list_to_string(lines[i:i+[x for x, val in enumerate(lines[i:i+20]) if '}' in val][0]]) for i, s in enumerate(lines) if ('RPKM' in s) and ('expression' in s) and ('}' not in lines[i+1])]
        
            locations_form = []
        
            for txt in locations:
                txt = txt.split()
                RPKM = [float(s[:s.find(')')]) for i, s in enumerate(txt) if ('RPKM' in txt[i-1])]
                description = [s for i, s in enumerate(txt) if ('text' not in s) and ('RPKM' not in s) and (')' not in s) and ('and' not in s) and ('tissues' not in s) and ('other' not in s) and (s.isdigit() == False)]
            
                locations_form.append((description,RPKM))
                
            print('\n',functions)
            print('\n',processes)
            print('\n',components)
            print('\n',locations_form)
            print(cluster)
            
        except:
            count = 0
            for ID in terms:
                count+=1
                handle = Entrez.efetch(db='gene', id=f'{ID} [mus musculus]', retmode = 'html')
                lines = handle.readlines()
                
                functions = [s.rstrip() for i,s in enumerate(lines) if ('enables' in lines[i-1]) and ('anchor' in s)]
                functions = np.unique(convert(functions), return_counts=True) 
                processes = [s.rstrip() for i,s in enumerate(lines) if ('involved' in lines[i-1]) and ('anchor' in s)]
                processes = np.unique(convert(processes), return_counts=True)
                components = [s.rstrip() for i,s in enumerate(lines) if ('located' in lines[i-1]) and ('anchor' in s)]
                components = np.unique(convert(components), return_counts=True)
            
                locations = [s.rstrip() for i, s in enumerate(lines) if ('RPKM' in s) and ('expression' in s) and ('}' in lines[i+1])]
                locations += [list_to_string(lines[i:i+[x for x, val in enumerate(lines[i:i+20]) if '}' in val][0]]) for i, s in enumerate(lines) if ('RPKM' in s) and ('expression' in s) and ('}' not in lines[i+1])]
            
                locations_form = []
                print(f"cluster {cluster}: {int((count/len(terms))*100)}%", end=' \r')
                
                for txt in locations:
                    txt = txt.split()
                    RPKM = [float(s[:s.find(')')]) for i, s in enumerate(txt) if ('RPKM' in txt[i-1])]
                    description = [s for i, s in enumerate(txt) if ('text' not in s) and ('RPKM' not in s) and ('and' not in s) and ('tissues' not in s) and ('other' not in s)]
                
                print(description)
            # print('\n',functions)
            # print('\n',processes)
            # print('\n',components)
            # print('\n',locations_form)
            # print(cluster)
            
        return locations_form

functions = get_location()
