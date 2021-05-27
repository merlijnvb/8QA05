import numpy as np

"""

   DEFINITION TO IMPORT TEXTFILES AND RETURN DICTIONARIES:
        
"""

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


"""

   K-MEANS CLUSTERING:
        
"""


class KMCA:
    def __init__(self, k=6, seed=2, data=dict(), normalize=True):
        '''
        Preconditions:  k (#nr of clusters), seed (#nr linked to the randomness), data (dictionary)
        Postconditions: save every datafield with the corresponding values:
                        --> self.k = #nr of clusters (k) .
                        --> self.seed = seed.
                        --> self.lib_[...] = empty dictionary that is filled later.
                            --> lib_data can be an empty dicionary or can be filled from the beginning. This choice is up to the user.
                        --> self.E_score = set to infinity, because if we don't cluster the data in the first iteration the E_score is not representative.
                        --> self.bool_normalized = boolean to check if the data is already normalized.
                        --> self.bool_want_normalized = boolean to check if user wants to normilize.
        '''
        
        self.k = k
        self.seed = seed
        self.lib_data = data
        self.lib_clustered = dict()
        self.lib_centroid = dict()
        self.E_score = float('inf')
        self.bool_normalized = False
        self.bool_want_normalized = normalize
    
    def normalize(self, data_unnorm):
        '''
        Preconditions:  data that is not yet normalized
        Postconditions: A dictionary:
                        --> Key = protein ID
                        --> Value = list of floats that are normalized by the use of the magnitude of the original vector to create a new vector with length 1
        Task of function: normalizing the data of every protein so it is a vector with length 1
        '''
        
        for index in data_unnorm:
            magnitude = np.sqrt(np.sum(np.square(np.array(data_unnorm[index]).astype(float)))) ## SET DATA TYPE TO FLOAT --> THEN CALCULATE THE LENGTH OF THE VECTOR BY APPLYING THE FORMULA IN THE CASUS
            self.lib_data[index] = np.divide(np.array(data_unnorm[index]).astype(float), magnitude) ## DIVIDE EACH VALUE IN THE CORDINATES BY THE LENGTH OF THE VECTOR
        
        self.bool_normalized = True
        
        return self.lib_data
    
    
    
    def cluster0(self):
        '''
        Postconditions: A dictionary:
                        --> Key = Cluster #nr
                        --> Value = list of protein ID's that are assigned to that cluster
        Task of function: assigning proteins to a random cluster
        '''
        
        # ERROR CHECHING --> WHEN FUNCTION IS CALLED BY USER AND NOT BY ALGORITHM ITSELF --> PREVENTING LARGER ERRORS
        if self.lib_data == dict(): # CHECKING IF LIB_DATA IS EMPTY
            raise ImportError("No input given. --> call function clustering(data) or __init__(data) or normalize(data)")
        
        np.random.seed(self.seed)
        
        labels = np.random.randint(0,self.k,len(self.lib_data))
        lib_labeled = dict()
        
        for i in range(len(self.lib_data)):
            lib_labeled[list(self.lib_data.keys())[i]] = labels[i]
            
        for k in range(self.k):
            self.lib_clustered[k] = np.array([i for i,j in lib_labeled.items() if j == k])
        
        return self.lib_clustered
    
    
    
    def centroid(self):
        '''
        Postconditions: A dictionary:
                        --> Key = Cluster #nr
                        --> Value = list containing the mean cordinate for every axis (list has the length of the #nr of dimensions the input has)
        Task of function: calculating the mean cordinate of every cluster (=centroid/centrum of every cluster)
        '''
        
        # ERROR CHECHING --> WHEN FUNCTION IS CALLED BY USER AND NOT BY ALGORITHM ITSELF --> PREVENTING LARGER ERRORS
        if self.lib_data == dict(): # CHECKING IF LIB_DATA IS EMPTY
            raise ImportError("No input given. --> clustering(data) or __init__(data) or normalize(data)")
            if self.lib_clustered == dict(): # CHECKING IF LIB_CLUSTERED IS EMPTY
                raise ImportError("No input given and no lib_clustered known. --> call function clustering(data) or __init__(data) or normalize(data) and call function cluster0()")
        if (self.lib_clustered == dict()) & (self.lib_data != dict()): # CHECKING IF LIB_CLUSTERED IS EMPTY BUT LIB_DATA IS GIVEN
            raise ImportError("No lib_clustered known. --> call function cluster0()")
            
        for k in self.lib_clustered: # CALCULATE FOR EVERY CLUSTER THE MEAN CORDINATE (=CENTROID)
            self.lib_centroid[k] = np.mean(np.array([self.lib_data[index] for index in self.lib_clustered[k]]), axis=0) 
            
        return self.lib_centroid
    
    
    
    def Escore(self):
        '''
        Postconditions: E_score => score how good the k-means clustering fit is
        Task of function: calculation the goodness of the k-means clustering application of certain centroids
        '''
        
        # ERROR CHECHING --> WHEN FUNCTION IS CALLED BY USER AND NOT BY ALGORITHM ITSELF --> PREVENTING LARGER ERRORS
        if self.lib_data == dict(): # CHECKING IF LIB_DATA IS EMPTY
            raise ImportError("No input given. --> call function clustering(data) or __init__(data) or normalize(data)")
            if self.lib_centroid == dict(): # CHECKING IF LIB_CENTROID IS EMPTY
                raise ImportError("No input given and no lib_centroid known. --> call function clustering(data) or __init__(data) or normalize(data) and call function centroid() or function clustering()")
        if (self.lib_centroid == dict()) & (self.lib_data != dict()): # CHECKING IF LIB_CENTROID IS EMPTY BUT LIB_DATA IS GIVEN
            raise ImportError("No lib_centroid known. --> call function centroid()")
            
        summation = 0
    
        for k in range(self.k): # CALCULATE FOR EVERY CLUSTER THE E_SCORE
            summation += np.sum([np.square(abs(np.subtract(self.lib_centroid[k], self.lib_data[dis]))) for dis in self.lib_clustered[k]])
            
        # CALCULATE THE MEAN E_SCORE FOR THE ENTIRE K-MEANS FIT
        E_score = summation / self.k
        
        return E_score
    
    
    def assign_cluster(self):
        '''
        Postconditions: A dictionary:
                        --> Key = Cluster #nr
                        --> Value = list of protein ID's that are assigned to the closest centroid
        Task of function: assigning proteins to the cluster/centroid which is closest to the protein
        '''
        
        # ERROR CHECHING --> WHEN FUNCTION IS CALLED BY USER AND NOT BY ALGORITHM ITSELF --> PREVENTING LARGER ERRORS
        if self.lib_data == dict(): # CHECKING IF LIB_DATA IS EMPTY
            raise ImportError("No input given. --> call function clustering(data) or __init__(data) or normalize(data)")
            if self.lib_centroid == dict(): # CHECKING IF LIB_CENTROIDS IS EMPTY
                raise ImportError("No input given and no lib_centroid known. --> call function clustering(data) or __init__(data) or normalize(data) and call function centroid() or function clustering()")
        if (self.lib_centroid == dict()) & (self.lib_data != dict()): # CHECKING IF LIB_CENTROID IS EMPTY BUT LIB_DATA IS GIVEN
            raise ImportError("No lib_centroid known. --> call function centroid()")
        
        self.lib_clustered = {k:[] for k in range(self.k)} # EMPTY THE CLUSTERED LIBRARY
        
        for index in self.lib_data: # LOOP OVER EVERY PROTEIN IN THE SELF.LIB_DATA
            # CALCULATE THE DISTANCE FROM THE PROTEIN TO EVERY CENTROID (CLUSTER) AND PICK THE CENTROID WHICH IS CLOSEST TO THE PROTEIN
            cluster = np.argmin([np.sqrt(np.sum(np.square(np.subtract(self.lib_centroid[k],self.lib_data[index])))) for k in range(self.k)])
            # APPEND THE CLUSTER LIST IN THE LIBRARY WITH THE CORRESPONDING PROTEIN ID
            self.lib_clustered[cluster].append(index)

        return self.lib_clustered
    
    
    
    def clustering(self, data=dict()):
        '''
        Postconditions: A dictionary:
                        --> Key = Cluster #nr
                        --> Value = list of protein ID's that are assigned to that cluster when the E_score (the goodness of the fit) is in its maximum
        Task of function: assigning proteins to the best cluster by updating the centroids a few iterations until the E_score is in its maximum
        '''
        
        # ERROR CHECHING --> WHEN FUNCTION IS CALLED BY USER AND NOT BY ALGORITHM ITSELF --> PREVENTING LARGER ERRORS
        if self.lib_data == dict(): # CHECKING IF LIB_DATA IS EMPTY
            if data == dict(): # CHECKING IF DATA IS EMPTY
                raise ImportError("No input given. --> clustering(data) or __init__(data) or normalize(data)")
            else:
                if self.bool_want_normalized: # CHECKING IF USER WANTS TO NORMALIZE
                    if self.bool_normalized == False: # CHECKING IF DATA ALREADY IS NORMALIZED
                        self.normalize(data)
                else:
                    self.lib_data = data
        if (self.bool_want_normalized) & (self.bool_normalized == False): # CHECKING IF LIB_DATA ALREADY IS GIVEN, BUT IS NOT NORMALIZED YET AND USER WANTS TO NORMALIZE (IS MOSTLY USED IF USER WANTS TO DIRECTLY OPTIMIZE THEIR CLUSTERING ALGORITHM)
            self.normalize(self.lib_data)
                 
        self.cluster0()
        self.centroid()
        self.E_score = self.Escore()
        
        optimized = False
        while optimized == False: # KEEP CALCULATING NEW CENTROIDS AND THEREFORE NEW CLUSTERS UNTIL THE E_SCORE IS AT ITS MAXIMUM
            self.lib_clustered = self.assign_cluster()
            self.lib_centroid = self.centroid()
            E_score_new = self.Escore()
            
            if self.E_score > E_score_new:
                self.E_score = E_score_new
            else:
                optimized = True
                
        return self.lib_clustered
        
    
    
    def optimize(self, min_seed=0, max_seed=30):
        '''
        Postconditions: A dictionary:
                        --> Key = Cluster #nr
                        --> Value = list of protein IDs that are assigned to that cluster when the E_score (the goodness of the fit) is in its maximum and when the starting centroids are at its best
        Task of function: calculating which seed gives the best results --> returning the clustered proteins for which the E_score is at its highest when comparing different starting centroids
        '''
        
        E_score_old = float('inf') # SET E_SCORE TO INFINITY SO THE 0TH ITERATION HAS THE LEAST FAVORABLE SCORE
        
        for seed in range(min_seed, max_seed+1): # LOOP OVER THE RANGE OF SEEDS THAT IS GIVEN
            self.seed = seed
            self.clustering()
            
            if self.E_score < E_score_old: # CHECK IF THE GIVEN SEED HAS A BETTER FIT THAN THE LAST BEST ONE, IFSO SAVE IT
                lib_best_clustered = self.lib_clustered
                best_seed = self.seed
                E_score_old = self.E_score
        
        # SAVE THE BEST RESULTS IN THE DATAFRAMES
        self.lib_clustered = lib_best_clustered
        self.seed = best_seed
        
        return self.lib_clustered
        
lib_data = file_to_lib('Data\Voorbeeld_clusterdata.txt')
lib_results = file_to_lib('Data\Voorbeeld_clusterresult.txt')
    
kmca = KMCA(data=lib_data)
#kmca.normalize(lib_data)
#res = kmca.clustering(lib_data)
best_res = kmca.optimize()

print(best_res, "Best seed:" ,kmca.seed)