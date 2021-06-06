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
    def __init__(self, k=6, seeds=[], data={}, normalize=True):
        '''
        Preconditions:  k (#nr of clusters), seed (#nr linked to the randomness), data (dictionary)
        Postconditions: save every datafield with the corresponding values:
                        --> self.k = #nr of clusters (k) .
                        --> self.seed = a list of all seeds to use for choosing random indexes for initial clustering.
                        --> self.bool_normalized = boolean to check if the data is already normalized.
                        --> self.bool_want_normalized = boolean to check if user wants to normilize.
                        --> self.lib_[...] = empty dictionary that is filled later.
                            --> lib_data can be an empty dicionary or can be filled from the beginning. This choice is up to the user.
                        --> self.E_score = set to infinity, because if we don't cluster the data in the first iteration the E_score is not representative.
                        --> self.Sil_score = set to 0, because if we don't cluster the data in the first iteration the Sil_score is not representative.
                        --> self.lib_Escores = empty dictionary that is later filled with different E-scores when parameters are changed.
                        --> self.lib_silscores = empty dictionary that is later filled with different Silhouette scores when parameters are changed.    
        '''
        self.k = k
        self.seeds = seeds
        if seeds != []:
            assert all([type(i)==int for i in seeds]),"The list 'seeds' may only contain integers as its values. --> __init__(seeds)"
        self.bool_normalized = False
        self.bool_want_normalized = normalize
        self.lib_data = data
        self.lib_clustered = {}
        self.lib_centroid = {}
        self.E_score = float('inf')
        self.Sil_score = -2
        self.lib_Escores = {}
        self.lib_silscores = {}

    
    
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
    
    
    
    def cluster0(self, seed=None):
        '''
        Preconditions: seed is an integer (the seed number to be used by the numpy.random module)
        Postconditions: A dictionary:
                        --> Key = Cluster #nr
                        --> Value = list of protein ID's that are assigned to that cluster
        Task of function: assigning proteins to a random cluster
        '''
        
        # ERROR CHECHING --> WHEN FUNCTION IS CALLED BY USER AND NOT BY ALGORITHM ITSELF --> PREVENTING LARGER ERRORS
        assert self.lib_data != {},"No data given. --> call function clustering(data) or __init__(data) or normalize(data)" # CHECKING IF LIB_DATA IS EMPTY
        
        np.random.seed(seed)
        
        randnr = np.random.randint(0,self.k,len(self.lib_data))
        labels = randnr#np.random.randint(0,self.k,len(self.lib_data))
        lib_labeled = {}
        
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
        
        # ERROR CHECKING --> WHEN FUNCTION IS CALLED BY USER AND NOT BY ALGORITHM ITSELF --> PREVENTING LARGER ERRORS
        if self.lib_data == {}: # CHECKING IF LIB_DATA IS EMPTY
            raise ImportError("No input given. --> clustering(data) or __init__(data) or normalize(data)")
            if self.lib_clustered == {}: # CHECKING IF LIB_CLUSTERED IS EMPTY
                raise ImportError("No input given and no lib_clustered known. --> call function clustering(data) or __init__(data) or normalize(data) and call function cluster0()")
        if (self.lib_clustered == {}) & (self.lib_data != {}): # CHECKING IF LIB_CLUSTERED IS EMPTY BUT LIB_DATA IS GIVEN
            raise ImportError("No lib_clustered known. --> call function cluster0()")
            
        for k in self.lib_clustered: # CALCULATE FOR EVERY CLUSTER THE MEAN CORDINATE (=CENTROID)
            self.lib_centroid[k] = np.mean(np.array([self.lib_data[index] for index in self.lib_clustered[k]]), axis=0) 
            
        return self.lib_centroid
    
    
    
    def Escore(self):
        '''
        Postconditions: E_score => score how good the k-means clustering fit is
        Task of function: calculation the goodness of the k-means clustering application of certain centroids
        '''
        
        # ERROR CHECKING --> WHEN FUNCTION IS CALLED BY USER AND NOT BY ALGORITHM ITSELF --> PREVENTING LARGER ERRORS
        if self.lib_data == {}: # CHECKING IF LIB_DATA IS EMPTY
            raise ImportError("No input given. --> call function clustering(data) or __init__(data) or normalize(data)")
            if self.lib_centroid == {}: # CHECKING IF LIB_CENTROID IS EMPTY
                raise ImportError("No input given and no lib_centroid known. --> call function clustering(data) or __init__(data) or normalize(data) and call function centroid() or function clustering()")
        if (self.lib_centroid == {}) & (self.lib_data != {}): # CHECKING IF LIB_CENTROID IS EMPTY BUT LIB_DATA IS GIVEN
            raise ImportError("No lib_centroid known. --> call function centroid()")
            
        summation = 0
    
        for k in range(self.k): # CALCULATE FOR EVERY CLUSTER THE E_SCORE
            summation += np.sum([np.square(abs(np.subtract(self.lib_centroid[k], self.lib_data[dis]))) for dis in self.lib_clustered[k]])
            
        # CALCULATE THE MEAN E_SCORE FOR THE ENTIRE K-MEANS FIT
        self.E_score = summation / self.k
        
        return self.E_score
    
    
    
    
    def silhouette_score(self):
        '''
        Postconditions: silhouette_score => score how good the k-means clustering fit is
        Task of function: calculating a score using a better scoring system with more indicative score values:
            --> A small internal dissimilarity value means it is well matched. Furthermore, a large external dissimilarity value means it is badly matched to its neighbouring cluster.
            
            --> Therefor if silhouette_score is closer to -1     ==>     datapoint would be better assigned to another cluster
            --> Therefor if silhouette_score is closer to  1     ==>     datapoint is appropriatly assigned to cluster
            --> If silhouette_score is close to 0     ==>     datapoint is on the border of two clusters.
        '''
        list_scores = []
        
        for index in self.lib_data:
            # GET THE CLUSTER #NR WHICH THE PROTEIN IS ASSIGNED TO
            cluster_nr = [cluster for cluster, indices in self.lib_clustered.items() if index in indices] # GET CLUSTER #NR THE PROTEIN IS ASSIGNED TO

            #CHECK IF THE PROTEIN IS ASSIGNED TO A CLUSTER
            if cluster_nr != []:
                cluster_nr = cluster_nr[0]
            
                if len(self.lib_clustered[cluster_nr]) > 1:
                    # CALCULATE THE MEAN DISTANCE OF THE PROTEIN TO ALL THE OTHER PROTEINS IN ITS CLUSTER (INTERNAL):
                    dissimilarity_internal = np.sum([np.linalg.norm(self.lib_data[index] - self.lib_data[index_in]) for index_in in self.lib_data]) / (len(self.lib_data)-1)
                
                    # CALCULATE THE MINIMAL MEAN DISTANCE OF THE PROTEIN TO ALL THE OTHER PROTEINS FROM THE OTHER CLUSTERS (EXTERNAL):
                    # --> MEAN DISTANCE IS GROUPED PER CLUSTER. FROM THIS LIST THE MINIMAL DISTANCE IS CALCULATED
                    dissimilarity_external = np.min([np.mean([np.linalg.norm(self.lib_data[index] - self.lib_data[index_ex]) for index_ex in self.lib_data if index_ex in self.lib_clustered[cluster]]) for cluster in self.lib_clustered if cluster != cluster_nr])
                
                    # CALCULATING THE SILHOUETTE SCORE FOR PROTEIN BY APPLYING THE FORMULA:
                    sil_score_i = (dissimilarity_external - dissimilarity_internal) / max(dissimilarity_external,dissimilarity_internal)
                
                if len(self.lib_clustered[cluster_nr]) == 1:
                    sil_score_i = 0
            
                list_scores.append(sil_score_i)
            
        self.Sil_score = np.mean(list_scores)
        
        if np.isnan(self.Sil_score):
            self.Sil_score = float('-inf')
        
        return self.Sil_score      
    
    
    def assign_cluster(self):
        '''
        Postconditions: A dictionary:
                        --> Key = Cluster #nr
                        --> Value = list of protein ID's that are assigned to the closest centroid
        Task of function: assigning proteins to the cluster/centroid which is closest to the protein
        '''
        
        # ERROR CHECKING --> WHEN FUNCTION IS CALLED BY USER AND NOT BY ALGORITHM ITSELF --> PREVENTING LARGER ERRORS
        if self.lib_data == {}: # CHECKING IF LIB_DATA IS EMPTY
            raise ImportError("No input given. --> call function clustering(data) or __init__(data) or normalize(data)")
            if self.lib_centroid == {}: # CHECKING IF LIB_CENTROIDS IS EMPTY
                raise ImportError("No input given and no lib_centroid known. --> call function clustering(data) or __init__(data) or normalize(data) and call function centroid() or function clustering()")
        if (self.lib_centroid == {}) & (self.lib_data != {}): # CHECKING IF LIB_CENTROID IS EMPTY BUT LIB_DATA IS GIVEN
            raise ImportError("No lib_centroid known. --> call function centroid()")
        
        self.lib_clustered = {k:[] for k in range(self.k)} # EMPTY THE CLUSTERED LIBRARY
        
        for index in self.lib_data: # LOOP OVER EVERY PROTEIN IN THE SELF.LIB_DATA
            # CALCULATE THE DISTANCE FROM THE PROTEIN TO EVERY CENTROID (CLUSTER) AND PICK THE CENTROID WHICH IS CLOSEST TO THE PROTEIN
            cluster = np.argmin([np.sqrt(np.sum(np.square(np.subtract(self.lib_centroid[k],self.lib_data[index])))) for k in range(self.k)])
            # APPEND THE CLUSTER LIST IN THE LIBRARY WITH THE CORRESPONDING PROTEIN ID
            self.lib_clustered[cluster].append(index)

        return self.lib_clustered
    
    
    
    def clustering(self, data={}, seed=None):
        '''
        Preconditions:  data: a dictionary:
                            --> Key = protein ID
                            --> Value = list of floats
                        seeds: a list of integers (the seed numbers to be used by the numpy.random module)
        Postconditions: A dictionary:
                        --> Key = cluster #nr
                        --> Value = list of protein ID's that are assigned to that cluster when the E_score (the goodness of the fit) is in its maximum
        Task of function: assigning proteins to the best cluster by updating the centroids a few iterations until the E_score is in its maximum
        '''
        
        # ERROR CHECKING --> WHEN FUNCTION IS CALLED BY USER AND NOT BY ALGORITHM ITSELF --> PREVENTING LARGER ERRORS
        assert type(data) == dict, "Input of 'data' is of a type other than dictionary, must be dictionary. --> clustering(data)"
        assert data != {} or self.lib_data != {}, "No input given or empty dictionary given as input for 'data'. --> optimize(data) or clustering(data) or __init__(data) or interval_data(data)"
        if data == {} and self.lib_data != {}:
            data = self.lib_data
        elif self.lib_data != data:
            self.lib_data = data
            
        if (self.bool_want_normalized) & (self.bool_normalized == False): # CHECKING IF LIB_DATA ALREADY IS GIVEN, BUT IS NOT NORMALIZED YET AND USER WANTS TO NORMALIZE (IS MOSTLY USED IF USER WANTS TO DIRECTLY OPTIMIZE THEIR CLUSTERING ALGORITHM)
            self.normalize(self.lib_data)
        
        self.cluster0(seed)
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
        
    
    def optimize(self, data={}, k_min=2, k_max=10, measure = 'Sil', seeds = []):
        '''
        Preconditions:  data: a dictionary:
                            --> Key = protein ID
                            --> Value = list of floats
                        k_min, k_max: range in which parameter k can be varied (lowest = k_min, highest = k_max)
                        measure: measure by which to evaluate how well the clustering fits (either 'E' for E-score or 'Sil' for silhouette score)
        Postconditions: lib_clustered: a dictionary:
                                    --> Key = Cluster #nr
                                    --> Value = list of protein IDs that are assigned to that cluster when the E_score (the goodness of the fit) is in its maximum
                        lib_Silscores/lib_Escores: a dictionary:
                                                --> Key = k
                                                --> Value = Sil_score/E_score
        Task of function: calculating which k gives the best results --> returning the clustered proteins for which the Sil_score is at its highest when comparing different k
        '''
        assert k_min < k_max, "k_min must be smaller than k_max. --> optimize(k_min,k_max) \n*If you only want to cluster the data for one value of k, call clustering()"
        assert k_min >= 1, 'k_min must be 1 or higher. --> optimize(k_min)'
        assert (measure == 'Sil' or measure == 'E'), "Measure must either be 'Sil' or 'E'. --> optimize(measure)"
        assert type(data) == dict, "Input of 'data' is of a type other than dictionary, must be dictionary. --> clustering(data)"
        assert data != {} or self.lib_data != {}, "No input given or empty dictionary given as input for 'data'. --> optimize(data) or clustering(data) or __init__(data) or interval_data(data)"
        
        if data == {} and self.lib_data != {}:
            data = self.lib_data
        elif self.lib_data != data:
            self.lib_data = data
            
        if seeds == []: 
            if self.seeds != []:
                seeds = self.seeds # IF THE LIST OF SEEDS IS NOT EXPLICITLY GIVEN TO THIS FUNCTION BUT ASSIGNED TO 'self.seeds' AT __init__, USE THE LIST self.seeds
            else:
                seeds = [None] # IF NO SEEDS ARE GIVEN TO THIS FUNCTION AND NONE WERE GIVEN AT self.seed, NO SEED SHOUL BE USED AT ALL (SO seeds MUST BE A LIST CONTAINING None AS ITS ONLY VALUE)
        assert all([type(i)==int for i in seeds]),"The list 'seeds' may only contain integers as its values. --> optimize(seeds)"
        
        for k in range(k_min, k_max+1): # LOOP OVER THE RANGE OF k'S THAT IS GIVEN
            self.k = k
            results = {} # THIS DICTIONARY WILL KEEP TRACK OF THE RESULTS FOUND FOR EVERY DIFFERENT SEED
            for seed in seeds: # THE CURRENT VALUE FOR k MUST BE TESTED USING EVERY SEED FROM THE LIST seeds
                print('k is ',k,'seed is ',seed)
                self.clustering(data)
                
                if measure == 'Sil':
                    try: # CALCULATE THE SILHOUETTE SCORE AND SAVE IT
                        self.silhouette_score()
                        results[seed] = self.Sil_score
                    except: # WHEN GIVEN NaN (WHEN NOT ALL CLUSTERS ARE FILLED OR ONLY 1 CLUSTER EXISTS) RETURN A VALUE THAT IS OUT OF THE SILHOUETTE RANGE
                        results[seed] = -2
                        
                else: # If not 'Sil', the measure is 'E' (see assert commands above)
                    self.Escore()
                    results[seed] = self.E_score
                print(results)  
            if measure == 'Sil':
                self.lib_silscores[k] = list(results.values())[np.argmax(list(results.values()))]
                print(list(results.keys()),list(results.values()),' of which ',self.lib_silscores[k],' is best')
            else:
                self.lib_Escores[k] = list(results.values())[np.argmin(list(results.values()))]
                print(list(results.keys()),list(results.values()),' of which ',self.lib_Escores[k],' is best')
                    
            print(f'{((k-k_min)/(k_max-k_min))*100}% ', end='\r') # PRINT AT EVERY ITERATION HOW FAR THE EVALUATION PROCES IS
        
        if measure == 'Sil':
            self.k =  list(self.lib_silscores.keys())[np.argmax(list(self.lib_silscores.values()))]
        else:
            self.k = list(self.lib_Escores.keys())[np.argmin(list(self.lib_Escores.values()))]
        
        self.clustering(data)
        self.Escore()
        self.silhouette_score()
        
        return self.lib_clustered



lib_data = file_to_lib('Data\Voorbeeld_clusterdata.txt')
lib_results = file_to_lib('Data\Voorbeeld_clusterresult.txt')
   
kmca = KMCA(data=lib_data,seeds=[0,1])
kmca_results = kmca.optimize(k_min=4,k_max=8,measure='E')
print(kmca.lib_silscores, kmca.lib_Escores)
print(kmca.k)
#kmca_scores = kmca.lib_silscores


"""
    PRINT TEXT FILES FORMATTED/SORTED LIKE THE INPUT FILE
    
"""


def return_txt_file(data, name, format_data=lib_data):
    lib_unfolded = dict()
    score = False
    
    for cluster in data:
        if type(data[cluster]) == int:
            lib_unfolded[data[cluster]] = cluster
        else:
            try:
                for ID in data[cluster]:
                    lib_unfolded[ID] = cluster
            except:
                score = True
                break
    
    file = open(f'{name}results.txt', 'w')
    indices_lost = list()
    
    if score == False:
        for INDEX in format_data:
            try:
                line = str(INDEX) + " " + str(lib_unfolded[INDEX]) + "\n"
                file.write(line)
            except:
                indices_lost.append(INDEX)
    else:
        for k in data:
            line = str(k) + " " + str(data[k]) + "\n"            
            file.write(line)
            
    file.close()
    
    if len(indices_lost) > 0:
        print(f'Function: {name} --> these indices are lost: {indices_lost}')
    
#return_txt_file(kmca_results, 'kmca_')
#return_txt_file(kmca_scores, 'kmca_sil_')