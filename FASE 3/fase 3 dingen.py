from Bio import Entrez
import numpy as np

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
lib_clust = file_to_lib('Data\clusterresultaten.txt')
lib_codes, columns = file_to_lib('Data/accessionnumbers.txt', True)

def get_location(code='ensemble'):
    location = [i for i, element in enumerate(columns[1:]) if code in element][0]
    
    for cluster in np.unique(list(lib_clust.values())):
        terms = [value[location] for index, value in lib_codes.items() if (index in lib_clust.keys()) and (lib_clust[index] == cluster) and (len(value) == 3)]
        terms_not_found = [index for index, value in lib_codes.items() if (index in lib_clust.keys()) and (lib_clust[index] == cluster) and (len(value) != 3)]
        
        Entrez.email = 'd.devetzis@student.tue.nl'        
        handle = Entrez.efetch(db='gene', id=terms, retmode = 'html')
        lines = handle.readlines()
        
        functions = [lines[i+1] for i,s in enumerate(lines) if ('enables' in s)  and ('anchor' in lines[i+1])]
        processes = [lines[i+1] for i,s in enumerate(lines) if ('involved' in s)  and ('anchor' in lines[i+1])]
        components = [lines[i+1] for i,s in enumerate(lines) if ('located' in s)  and ('anchor' in lines[i+1])]
        locations = [lines[i+1:i+2] for i, s in enumerate(lines) if ('label "Text Summary"' in s)]
        
get_location()