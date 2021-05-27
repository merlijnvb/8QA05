import numpy as np

def file_to_lib(FileName):
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
        line = lines.split()
        if len(line[1:]) == 1:
            lib[int(line[0])] = int(line[1])
        else:
            lib[int(line[0])] = np.array(line[1:]).astype(float)
            
    return lib

class GCA:
    def __init__(self, scope=1, subspaces=9, tresh=3, data=dict()):
        self.scope = scope
        self.subspaces = subspaces 
        self.tresh = tresh 
        self.lib_data = data
        self.lib_axis = dict()
        self.lib_interval = dict()
        self.lib_cells = dict()
        
    def interval_data(self, data=dict()):
        if self.lib_data == dict():
            self.lib_data = data
        
        
        
        for axis in range(len(list(self.lib_data.values())[0])): 
            x = list()
            for index in self.lib_data: 
                x.append(float(self.lib_data[index][axis]))  
        
            self.lib_axis[axis] = np.array(x) 
       
        for axis in self.lib_axis:
            minimum = np.min(self.lib_axis[axis])
            maximum = np.max(self.lib_axis[axis])

            self.lib_interval[axis] = np.mgrid[minimum:maximum:complex(0,self.subspaces+1)]
        
    def grid_data(self):
        for axis in self.lib_interval: # LOOP OVER EVERY AXIS
            for interval in range(1,len(self.lib_interval[axis])): # LOOP OVER EVERY SUBSPACE (SPACE BETWEEN 2 INTERVALS)
                # CONDITIONS ARE APPLIED TO THE VALUES OF AN AXIS
                # --> CONDITIONS ARE THAT THE VALUES NEED TO BE BETWEEN 2 INTERVALS SET BY THE INTERVAL_DATA FUNCTION
                # --> LIB_CORDS => A LIST WITH AS VALUES THE SUBSPACE-ID 
                self.lib_axis[axis][(self.lib_axis[axis] >= self.lib_interval[axis][interval-1]) & (self.lib_axis[axis] <= self.lib_interval[axis][interval]) & (type(self.lib_axis[axis]) != str)] = str(interval)
            
            # APPEND THE GRID LIST WITH THE CALCULATED SUBSPACE LISTS PER AXIS
            print( in self.lib_axis[axis])
            #grid.append(self.lib_axis[axis])
                
        #for i in  self.lib_axis:
             #for j in self.lib_axis[i]:
                 #print(type(j))
        


lib_data = file_to_lib('Data\Voorbeeld_clusterdata.txt')
lib_results = file_to_lib('Data\Voorbeeld_clusterresult.txt')
   
gca = GCA()
gca.interval_data(lib_data)
gca.grid_data()