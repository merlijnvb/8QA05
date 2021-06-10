"""
Course code: 8QA05
group number: 9
members: Merlijn van Benthem, Dimo Devetzis en Sam Reijs
"""


import numpy as np


"""
   FUNCTION TO IMPORT TEXTFILES AND RETURN DICTIONARIES
        
"""
def file_to_dict(FileName):
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
        
    data_dict = {}
    
    for lines in opened_data:
        line = lines.split()
        if len(line[1:]) == 1:
            data_dict[int(line[0])] = int(line[1])
        else:
            data_dict[int(line[0])] = np.array(line[1:]).astype(float)
            
    return data_dict


"""
    K-MEANS CLUSTERING ALGORITHM (KMCA)
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
                        --> self.dict_[...] = empty dictionary that is filled later.
                            --> dict_data can be an empty dicionary or can be filled from the beginning. This choice is up to the user.
                        --> self.E_score = set to infinity, because if we don't cluster the data in the first iteration the E_score is not representative.
                        --> self.Sil_score = set to 0, because if we don't cluster the data in the first iteration the Sil_score is not representative.
                        --> self.dict_Escores = empty dictionary that is later filled with different E-scores when parameters are changed.
                        --> self.dict_silscores = empty dictionary that is later filled with different Silhouette scores when parameters are changed.    
        '''
        self.k = k
        self.seeds = seeds
        if seeds != []:
            assert all([type(i)==int for i in seeds]),"The list 'seeds' may only contain integers as its values. --> __init__(seeds)"
        self.bool_normalized = False
        self.bool_want_normalized = normalize
        self.dict_data = data
        self.dict_clustered = {}
        self.dict_centroid = {}
        self.E_score = float('inf')
        self.Sil_score = -2
        self.dict_Escores = {}
        self.dict_silscores = {}

    
    
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
            self.dict_data[index] = np.divide(np.array(data_unnorm[index]).astype(float), magnitude) ## DIVIDE EACH VALUE IN THE CORDINATES BY THE LENGTH OF THE VECTOR
        
        self.bool_normalized = True
        
        return self.dict_data
    
    
    
    def cluster0(self, seed=None):
        '''
        Preconditions: seed is an integer (the seed number to be used by the numpy.random module)
        Postconditions: A dictionary:
                        --> Key = Cluster #nr
                        --> Value = list of protein ID's that are assigned to that cluster
        Task of function: assigning proteins to a random cluster
        '''
        
        # ERROR CHECHING --> WHEN FUNCTION IS CALLED BY USER AND NOT BY ALGORITHM ITSELF --> PREVENTING LARGER ERRORS
        assert self.dict_data != {},"No data given. --> call function clustering(data) or __init__(data) or normalize(data)" # CHECKING IF dict_data IS EMPTY
        
        np.random.seed(seed)
        
        randnr = np.random.randint(0,self.k,len(self.dict_data))
        labels = randnr#np.random.randint(0,self.k,len(self.dict_data))
        dict_labeled = {}
        
        for i in range(len(self.dict_data)):
            dict_labeled[list(self.dict_data.keys())[i]] = labels[i]
            
        for k in range(self.k):
            self.dict_clustered[k] = np.array([i for i,j in dict_labeled.items() if j == k])
        
        return self.dict_clustered
    
    
    
    def centroid(self):
        '''
        Postconditions: A dictionary:
                        --> Key = Cluster #nr
                        --> Value = list containing the mean cordinate for every axis (list has the length of the #nr of dimensions the input has)
        Task of function: calculating the mean cordinate of every cluster (=centroid/centrum of every cluster)
        '''
        
        # ERROR CHECKING --> WHEN FUNCTION IS CALLED BY USER AND NOT BY ALGORITHM ITSELF --> PREVENTING LARGER ERRORS
        if self.dict_data == {}: # CHECKING IF dict_data IS EMPTY
            raise ImportError("No input given. --> clustering(data) or __init__(data) or normalize(data)")
            if self.dict_clustered == {}: # CHECKING IF dict_clustered IS EMPTY
                raise ImportError("No input given and no dict_clustered known. --> call function clustering(data) or __init__(data) or normalize(data) and call function cluster0()")
        if (self.dict_clustered == {}) & (self.dict_data != {}): # CHECKING IF dict_clustered IS EMPTY BUT dict_data IS GIVEN
            raise ImportError("No dict_clustered known. --> call function cluster0()")
            
        for k in self.dict_clustered: # CALCULATE FOR EVERY CLUSTER THE MEAN CORDINATE (=CENTROID)
            self.dict_centroid[k] = np.mean(np.array([self.dict_data[index] for index in self.dict_clustered[k]]), axis=0) 
            
        return self.dict_centroid
    
    
    
    def Escore(self):
        '''
        Postconditions: E_score => score how good the k-means clustering fit is
        Task of function: calculation the goodness of the k-means clustering application of certain centroids
        '''
        
        # ERROR CHECKING --> WHEN FUNCTION IS CALLED BY USER AND NOT BY ALGORITHM ITSELF --> PREVENTING LARGER ERRORS
        if self.dict_data == {}: # CHECKING IF dict_data IS EMPTY
            raise ImportError("No input given. --> call function clustering(data) or __init__(data) or normalize(data)")
            if self.dict_centroid == {}: # CHECKING IF dict_centroid IS EMPTY
                raise ImportError("No input given and no dict_centroid known. --> call function clustering(data) or __init__(data) or normalize(data) and call function centroid() or function clustering()")
        if (self.dict_centroid == {}) & (self.dict_data != {}): # CHECKING IF dict_centroid IS EMPTY BUT dict_data IS GIVEN
            raise ImportError("No dict_centroid known. --> call function centroid()")
            
        summation = 0
    
        for k in range(self.k): # CALCULATE FOR EVERY CLUSTER THE E_SCORE
            summation += np.sum([np.square(abs(np.subtract(self.dict_centroid[k], self.dict_data[dis]))) for dis in self.dict_clustered[k]])
            
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
        
        for index in self.dict_data:
            # GET THE CLUSTER #NR WHICH THE PROTEIN IS ASSIGNED TO
            cluster_nr = [cluster for cluster, indices in self.dict_clustered.items() if index in indices] # GET CLUSTER #NR THE PROTEIN IS ASSIGNED TO

            #CHECK IF THE PROTEIN IS ASSIGNED TO A CLUSTER
            if cluster_nr != []:
                cluster_nr = cluster_nr[0]
            
                if len(self.dict_clustered[cluster_nr]) > 1:
                    # CALCULATE THE MEAN DISTANCE OF THE PROTEIN TO ALL THE OTHER PROTEINS IN ITS CLUSTER (INTERNAL):
                    dissimilarity_internal = np.sum([np.linalg.norm(self.dict_data[index] - self.dict_data[index_in]) for index_in in self.dict_data]) / (len(self.dict_data)-1)
                
                    # CALCULATE THE MINIMAL MEAN DISTANCE OF THE PROTEIN TO ALL THE OTHER PROTEINS FROM THE OTHER CLUSTERS (EXTERNAL):
                    # --> MEAN DISTANCE IS GROUPED PER CLUSTER. FROM THIS LIST THE MINIMAL DISTANCE IS CALCULATED
                    dissimilarity_external = np.min([np.mean([np.linalg.norm(self.dict_data[index] - self.dict_data[index_ex]) for index_ex in self.dict_data if index_ex in self.dict_clustered[cluster]]) for cluster in self.dict_clustered if cluster != cluster_nr])
                
                    # CALCULATING THE SILHOUETTE SCORE FOR PROTEIN BY APPLYING THE FORMULA:
                    sil_score_i = (dissimilarity_external - dissimilarity_internal) / max(dissimilarity_external,dissimilarity_internal)
                
                if len(self.dict_clustered[cluster_nr]) == 1:
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
        if self.dict_data == {}: # CHECKING IF dict_data IS EMPTY
            raise ImportError("No input given. --> call function clustering(data) or __init__(data) or normalize(data)")
            if self.dict_centroid == {}: # CHECKING IF dict_centroidS IS EMPTY
                raise ImportError("No input given and no dict_centroid known. --> call function clustering(data) or __init__(data) or normalize(data) and call function centroid() or function clustering()")
        if (self.dict_centroid == {}) & (self.dict_data != {}): # CHECKING IF dict_centroid IS EMPTY BUT dict_data IS GIVEN
            raise ImportError("No dict_centroid known. --> call function centroid()")
        
        self.dict_clustered = {k:[] for k in range(self.k)} # EMPTY THE CLUSTERED DICTIONARY
        
        for index in self.dict_data: # LOOP OVER EVERY PROTEIN IN THE SELF.dict_data
            # CALCULATE THE DISTANCE FROM THE PROTEIN TO EVERY CENTROID (CLUSTER) AND PICK THE CENTROID WHICH IS CLOSEST TO THE PROTEIN
            cluster = np.argmin([np.sqrt(np.sum(np.square(np.subtract(self.dict_centroid[k],self.dict_data[index])))) for k in range(self.k)])
            # APPEND THE CLUSTER LIST IN THE DICTIONARY WITH THE CORRESPONDING PROTEIN ID
            self.dict_clustered[cluster].append(index)

        return self.dict_clustered
    
    
    
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
        assert data != {} or self.dict_data != {}, "No input given or empty dictionary given as input for 'data'. --> optimize(data) or clustering(data) or __init__(data) or interval_data(data)"
        if data == {} and self.dict_data != {}:
            data = self.dict_data
        elif self.dict_data != data:
            self.dict_data = data
            
        if (self.bool_want_normalized) & (self.bool_normalized == False): # CHECKING IF dict_data ALREADY IS GIVEN, BUT IS NOT NORMALIZED YET AND USER WANTS TO NORMALIZE (IS MOSTLY USED IF USER WANTS TO DIRECTLY OPTIMIZE THEIR CLUSTERING ALGORITHM)
            self.normalize(self.dict_data)
        
        self.cluster0(seed)
        self.centroid()
        self.E_score = self.Escore()
        
        optimized = False
        while optimized == False: # KEEP CALCULATING NEW CENTROIDS AND THEREFORE NEW CLUSTERS UNTIL THE E_SCORE IS AT ITS MAXIMUM
            self.dict_clustered = self.assign_cluster()
            self.dict_centroid = self.centroid()
            E_score_new = self.Escore()
                
            if self.E_score > E_score_new:
                self.E_score = E_score_new
            else:
                optimized = True
            
        return self.dict_clustered
        
    
    def optimize(self, data={}, k_min=2, k_max=10, measure = 'Sil', seeds = []):
        '''
        Preconditions:  data: a dictionary:
                            --> Key = protein ID
                            --> Value = list of floats
                        k_min, k_max: range in which parameter k can be varied (lowest = k_min, highest = k_max)
                        measure: measure by which to evaluate how well the clustering fits (either 'E' for E-score or 'Sil' for silhouette score)
        Postconditions: dict_clustered: a dictionary:
                                    --> Key = Cluster #nr
                                    --> Value = list of protein IDs that are assigned to that cluster when the E_score (the goodness of the fit) is in its maximum
                        dict_silscores/dict_Escores: a dictionary:
                                                --> Key = k
                                                --> Value = Sil_score/E_score
        Task of function: calculating which k gives the best results --> returning the clustered proteins for which the Sil_score is at its highest when comparing different k
        '''
        assert k_min <= k_max, "k_min must be smaller than or equal to k_max. --> optimize(k_min,k_max)"
        assert k_min >= 1, 'k_min must be 1 or higher. --> optimize(k_min)'
        assert (measure == 'Sil' or measure == 'E'), "Measure must either be 'Sil' or 'E'. --> optimize(measure)"
        assert type(data) == dict, "Input of 'data' is of a type other than dictionary, must be dictionary. --> clustering(data)"
        assert data != {} or self.dict_data != {}, "No input given or empty dictionary given as input for 'data'. --> optimize(data) or clustering(data) or __init__(data) or interval_data(data)"
        
        if data == {} and self.dict_data != {}:
            data = self.dict_data
        elif self.dict_data != data:
            self.dict_data = data
            
        if seeds == []: 
            if self.seeds != []:
                seeds = self.seeds # IF THE LIST OF SEEDS IS NOT EXPLICITLY GIVEN TO THIS FUNCTION BUT ASSIGNED TO 'self.seeds' AT __init__, USE THE LIST self.seeds
            else:
                seeds = [None] # IF NO SEEDS ARE GIVEN TO THIS FUNCTION AND NONE WERE GIVEN AT self.seed, NO SEED SHOUL BE USED AT ALL (SO seeds MUST BE A LIST CONTAINING None AS ITS ONLY VALUE)
        assert all([type(i)==int for i in seeds]),"The list 'seeds' may only contain integers as its values. --> optimize(seeds)"
        print('\nFinding the optimal number of clusters\nPROGRESS:')
        for k in range(k_min, k_max+1): # LOOP OVER THE RANGE OF k'S THAT IS GIVEN
            self.k = k
            results = {} # THIS DICTIONARY WILL KEEP TRACK OF THE RESULTS FOUND FOR EVERY DIFFERENT SEED
            for seed in seeds: # THE CURRENT VALUE FOR k MUST BE TESTED USING EVERY SEED FROM THE LIST seeds
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

            if measure == 'Sil':
                self.dict_silscores[k] = list(results.values())[np.argmax(list(results.values()))]
            else:
                self.dict_Escores[k] = list(results.values())[np.argmin(list(results.values()))]
                    
            if k_min != k_max: # IF k_min IS EQUAL TO k_max, ONLY ONE VALUE FOR k IS TESTED. IN THAT CASE, THE PROGRESS BAR IS NOT NECESSARY.
                print(f'{(((k-k_min)/(k_max-k_min))*100):.2f}%') # PRINT AT EVERY ITERATION HOW FAR THE EVALUATION PROCES IS
                
        if measure == 'Sil':
            self.k = list(self.dict_silscores.keys())[np.argmax(list(self.dict_silscores.values()))]
            self.Sil_score = self.dict_silscores[self.k]
            self.clustering()
            self.Escore()
        else:
            self.k = list(self.dict_Escores.keys())[np.argmin(list(self.dict_Escores.values()))]
            self.E_score = self.dict_Escores[self.k]
            self.clustering()
            self.silhouette_score()
        
        return self.dict_clustered



"""
    GRID-BASED CLUSTERING ALGORITHM (GBCA)
"""
class GBCA:
    def __init__(self, scope=1, subspaces=9, thres=3, data={}):
        '''
        Preconditions:  scope: range in which cells are labeled as neigbours
                        subspaces: #nr of subspaces there are for every axis
                        thres: threshold for the amount of data that a cluster needs to contain to be labeled as a cluster
        Postconditions: save every datafield with the corresponding values:
                        --> self.scope = range of neigbours (scope).
                        --> self.subspaces = #nr of subpaces (subspaces).
                        --> self.thres = #nr of datapoints within group of cells (thres).
                        --> self.dict_[...] = empty dictionary that is filled later.
                            --> dict_data can be an empty dicionary or can be filled from the beginning. This choice is up to the user.
                        --> self.indices_not_clustered = empty array that is filled later.
                        --> self.E_score = set to infinity, because if we don't cluster the data in the first iteration the E_score is not representative.
                        --> self.Sil_score = set to 0, because if we don't cluster the data in the first iteration the Sil_score is not representative.
                        --> self.dict_silscores = empty dictionary that is later filled with different Silhouette scores when parameters are changed.
                        --> self.dict_Escores = empty dictionary that is later filled with different E-scores when parameters are changed.
        '''
        self.scope = scope
        self.subspaces = subspaces 
        self.thres = thres 
        self.dict_data = data
        self.dict_axis = {}
        self.dict_interval = {}
        self.dict_grid = {}
        self.dict_clustered = {}    
        self.indices_not_clustered = np.array([])
        self.E_score = float('inf')
        self.Sil_score = -2
        self.dict_Escores = {}
        self.dict_silscores = {}
            
        
    def interval_data(self, data={}):
        '''
        Precondition: A dictionary:
                        --> Key = protein ID
                        --> Value = list of floats
        Postconditions: A dictionary:
                        --> Key = Axis
                        --> Value = list of floats which represent the intervals that make the boundaries of the subspaces
        Task of function: calculate the boundries of the supspaces for every axis
        '''
        assert type(data) == dict, "Input of 'data' is of a type other than dictionary, must be dictionary. --> clustering(data)"
        assert data != {} or self.dict_data != {}, "No input given or empty dictionary given as input for 'data'. --> optimize(data) or clustering(data) or __init__(data) or interval_data(data)"
        if data == {} and self.dict_data != {}:
            data = self.dict_data
        elif self.dict_data != data:
            self.dict_data = data
        
        # CONVERT THE DATA IN dict_data TO A DICTIONARY WITH THE KEYS BEING THE AXES AND THE VALUES BEING THE CORRESPONDING VALUES OF EVERY AXIS FOR EVERY PROTEIN
        for axis in range(len(list(self.dict_data.values())[0])): 
            x = []
            for index in self.dict_data: 
                x.append(float(self.dict_data[index][axis]))  
        
            self.dict_axis[axis] = np.array(x) 
       
        # CALCULATE THE INTERVALS PER AXIS USING THE MINIMUM AND MAXIMUM (RANGE IN WHICH ALL DATA POINTS LAY)
        for axis in self.dict_axis:
            minimum = np.min(self.dict_axis[axis])
            maximum = np.max(self.dict_axis[axis])

            self.dict_interval[axis] = np.mgrid[minimum:maximum:complex(0,self.subspaces+1)] # MAKE 1-D GRID => ((abs(min)+max/subspaces) * n_step) + min
        
        
        
    def grid_data(self):
        '''
        Postconditions: A dictionary:
                        --> Key = protein ID
                        --> Value = list of integers which represent the subspaces of the axes (the list represents the cell ID)
        Task of function: calculate the unique cell ID for every protein
        '''
        
        grid = {}
        
        for axis in self.dict_interval: # LOOP OVER EVERY AXIS
            dict_axis_subspace = self.dict_axis[axis].copy() # COPY THE AXIS LIST TO CONVERT THE DATAPOINTS WHICH FALLS BETWEEN INTERVALS TO SUBSPACE IDs
            
            for interval in range(1,len(self.dict_interval[axis])): # LOOP OVER EVERY INTERVAL EXCEPT THE FIRST ONE
                # APPLY THE INTERVAL CONDITIONS TO THE CORRESPONDING AXIS AND CONVERT THE DATAPOINTS WHICH RETURN TRUE TO THE CONDITIONS TO THE SUBSPACE ID
                dict_axis_subspace[(self.dict_axis[axis] >= self.dict_interval[axis][interval-1]) & (self.dict_axis[axis] <= self.dict_interval[axis][interval])] = interval
                         
            # SAVE A DICTIONARY WITH AS KEY THE AXIS AND AS VALUE A LIST OF THE ASSIGNED SUBSPACE IDs FOR EVERY PROTEIN
            grid[axis] = dict_axis_subspace.astype(int)
            
        for i in range(len(self.dict_data)): # LOOP OVER EVERY PROTEIN
            # CONVERT THE GRID DICTIONARY TO A DICTIONARY WITH THE KEYS SET TO THE PROTEIN ID AND THE VALUE SET TO A LIST CONTAINING ALL THE SUBSPACE IDs FOR EVERY AXIS ==> CELL ID
            index = list(self.dict_data.keys())[i]
            self.dict_grid[index] = np.zeros((len(self.dict_interval),),dtype=int)
            
            for axis in grid:
                self.dict_grid[index][axis] = grid[axis][i]        
                
    
    def assign_neigbours(self, index0, cluster):
        '''
        Precondition: 
                    index0: protein ID of the protein for which the neigbours need to be found
                    cluster: list to which the newfound neigbours are assigned
                    options: list of the ID's of all proteins that are not yet clustered (and thus possible neighbours to index0)
        Postconditions: the protein IDs of all the neigbouring proteins to the starting protein are added to the list 'cluster'
        Task of function: find all the neigbours of starting protein x and assign it to the same group, then repeat for those neighbours
        '''
        neighbours = []
    
        for index in self.indices_not_clustered: # LOOP OVER THE NOT YET ASSIGNED PROTEINS
            # CHECK FOR EVERY, NOT YET ASSIGNED, PROTEIN IF THEY ARE IN RANGE TO BE CALLED A NEIGBOUR
            # --> IF TRUE: GROUP IT WITH THE ALREADY EXISTING GROUP OR GROUP IT WITH THE ORIGINAL PROTEIN
            if all([-self.scope <= axis_distance <= self.scope for axis_distance in np.subtract(self.dict_grid[index], self.dict_grid[index0])]):
                neighbours.append(index)
                cluster.append(index)
                self.indices_not_clustered = np.delete(self.indices_not_clustered, np.where(self.indices_not_clustered == index))          
        
        # IF THE LIST CONTAINING NEIGBOURS CONSIST MORE THAN ONLY ITS OWN ID THEN IT CAN CALL ON ITS SELF TO FIND ALL THE OTHER NEIGBOURS UNTIL THERE ARE NONE TO BE FOUND
        if len(neighbours) > 1:
             for neighbour in neighbours:
                 self.assign_neigbours(neighbour,cluster)
        
        
    
    def clustering(self, data={}):
        '''
        Precondition: a dictionary:
                        --> Key = protein ID
                        --> Value = list of floats
        Postconditions: a dictionary:
                        --> Key = cluster #nr
                        --> Value = list of protein ID's that are assigned to that cluster and are neigbours from eachother
        Task of function: cluster all the neigbours of every protein not yet assigned to a cluster 
        '''
        assert type(data) == dict, "Input of 'data' is of a type other than dictionary, must be dictionary. --> clustering(data)"
        assert data != {} or self.dict_data != {}, "No input given or empty dictionary given as input for 'data'. --> optimize(data) or clustering(data) or __init__(data) or interval_data(data)"
        if data == {} and self.dict_data != {}:
            data = self.dict_data
        elif self.dict_data != data:
            self.dict_data = data
        
        self.indices_not_clustered = np.array(list(self.dict_data.keys()))
        self.interval_data()
        self.grid_data()
        
        cluster_nr = 0
        anticluster_nr = 0
        dict_anticlustered = {}
        
        for index in self.indices_not_clustered: # LOOP OVER EVERY PROTEIN THAT IS NOT YET ASSIGNED TO ANY CLUSTER OR ANTICLUSTER
            if index in self.indices_not_clustered: 
                cluster = []
                self.assign_neigbours(index,cluster)
                
                # CHECK IF THE CLUSTER IS BIGGER THEN THE SET THRESHOLD
                # --> IF TRUE: ADD THE CLUSTER TO THE DICTIONARY CONTAINING THE CLUSTERS
                # --> IF FALSE: ADD THE CLUSTER TO THE DICTIONARY CONTAINING THE NON-CLUSTERS IN THE DICTIONARY CONTAINING THE CLUSTERS
                if len(cluster) > self.thres:              
                    self.dict_clustered[cluster_nr] = cluster
                    cluster_nr += 1
                else:
                    dict_anticlustered[anticluster_nr] = cluster
                    anticluster_nr += 1
            
        # APPEND THE CLUSTER DICTIONARY WITH THE DICTIONARY CONTAINING THE NON-CLUSTERS             
        self.dict_clustered[-1] = dict_anticlustered
        
        return self.dict_clustered
    
    
    
    def Escore(self):
        '''
        Postconditions: E_score => score how good the grid-based clustering fit is
        Task of function: calculation the goodness of the k-means clustering application of certain centroids
        '''
        
        summation = 0
        dict_centroid = {}
        
        # CALCULATE THE CENTROIDS OF EVERY CLUSTER
        for k in self.dict_clustered: # LOOP OVER THE CLUSTERS
            if k >= 0: # ONLY CALCULATE THE CENTROIDS FOR THE CLUSTERS NOT FOR THE NOT ASSIGNED DATAPOINTS
                coordinate = []
                for axis in self.dict_axis: # LOOP OVER THE AXES
                    mean_x = 0
                    
                    for protein in [value for index, value in self.dict_data.items() if index in self.dict_clustered[k]]: # CALCULATE THE MEAN VALUE OF THE AXIS
                        mean_x += protein[axis]
                    
                    coordinate.append(mean_x / len(self.dict_clustered[k]))
                    dict_centroid[k] = coordinate
        
        for k in dict_centroid: # CALCULATE FOR EVERY CLUSTER THE E_SCORE
            summation += np.sum([np.square(abs(np.subtract(dict_centroid[k], self.dict_data[dis]))) for dis in self.dict_clustered[k]])
            
        # CALCULATE THE MEAN E_SCORE FOR THE ENTIRE GRID CLUSTERING FIT
        self.E_score = summation / len(dict_centroid)
        
        return self.E_score
    
    
    
    def silhouette_score(self):
        '''
        Postconditions: silhouette_score => score how good the grid-based clustering fit is
        Task of function: calculating a score using a better scoring system with more indicative score values:
            --> A small internal dissimilarity value means it is well matched. Furthermore, a large external dissimilarity value means it is badly matched to its neighbouring cluster.
            
            --> Therefor if silhouette_score is closer to -1     ==>     datapoint would be better assigned to another cluster
            --> Therefor if silhouette_score is closer to  1     ==>     datapoint is appropriatly assigned to cluster
            --> If silhouette_score is close to 0     ==>     datapoint is on the border of two clusters.
        '''
        list_scores = []
        
        for index in self.dict_data:
            # GET THE CLUSTER #NR WHICH THE PROTEIN IS ASSIGNED TO
            cluster_nr = [cluster for cluster, indices in self.dict_clustered.items() if (index in indices) and (cluster != -1)] # GET CLUSTER #NR THE PROTEIN IS ASSIGNED TO
            
            #CHECK IF THE PROTEIN IS ASSIGNED TO A CLUSTER
            if cluster_nr != []:
                cluster_nr = cluster_nr[0]
            
                if len(self.dict_clustered[cluster_nr]) > 1:
                    # CALCULATE THE MEAN DISTANCE OF THE PROTEIN TO ALL THE OTHER PROTEINS IN ITS CLUSTER (INTERNAL):
                    dissimilarity_internal = np.sum([np.linalg.norm(self.dict_data[index] - self.dict_data[index_in]) for index_in in self.dict_data]) / (len(self.dict_data)-1)
                    
                    # CALCULATE THE MINIMAL MEAN DISTANCE OF THE PROTEIN TO ALL THE OTHER PROTEINS FROM THE OTHER CLUSTERS (EXTERNAL):
                    # --> MEAN DISTANCE IS GROUPED PER CLUSTER. FROM THIS LIST THE MINIMAL DISTANCE IS CALCULATED
                    dissimilarity_external = np.min([np.mean([np.linalg.norm(self.dict_data[index] - self.dict_data[index_ex]) for index_ex in self.dict_data if index_ex in self.dict_clustered[cluster]]) for cluster in self.dict_clustered if (cluster != cluster_nr) & (cluster != -1)])
                    
                    # CALCULATING THE SILHOUETTE SCORE FOR PROTEIN BY APPLYING THE FORMULA:
                    sil_score_i = (dissimilarity_external - dissimilarity_internal) / max(dissimilarity_external,dissimilarity_internal)
                    
                if len(self.dict_clustered[cluster_nr]) == 1:
                    sil_score_i = 0
                
                list_scores.append(sil_score_i)
            
        self.Sil_score = np.mean(list_scores)
        
        if np.isnan(self.Sil_score):
            self.Sil_score = float('-inf')
                
        return self.Sil_score
    
    
    def optimize(self, data={}, subspace_min=10, subspace_max=20, measure = 'Sil'):
        '''
        Preconditions:  data: a dictionary:
                            --> Key = protein ID
                            --> Value = list of floats
                        subspace_min, subspace_max: range in which parameter subspaces can be varied (lowest = subspace_min, highest = subspace_max)
                        measure: measure by which to evaluate how well the clustering fits (either 'E' for E-score or 'Sil' for silhouette score)
        Postconditions: dict_clustered: a dictionary:
                                    --> Key = Cluster #nr
                                    --> Value = list of protein IDs that are assigned to that cluster when the E_score (the goodness of the fit) is in its maximum
                        dict_silscores/dict_Escores: a dictionary:
                                                --> Key = subspace
                                                --> Value = Sil_score/E_score
        Task of function: calculating which subspace gives the best results --> returning the clustered proteins for which the Sil_score is at its highest when comparing different subspaces
        '''
        assert subspace_min < subspace_max, "subspace_min must be smaller than subspace_max. --> optimize(subspace_min,subspace_max)"
        assert subspace_min >= 1, "subspace_min must be 1 or higher. --> optimize(subspace_min)"
        assert (measure == 'Sil' or measure == 'E'), "Measure must either be 'Sil' or 'E'. --> optimize(measure)"
        assert type(data) == dict, "Input of 'data' is of a type other than dictionary, must be dictionary. --> clustering(data)"
        assert data != {} or self.dict_data != {}, "No input given or empty dictionary given as input for 'data'. --> optimize(data) or clustering(data) or __init__(data) or interval_data(data)"
        if data == {} and self.dict_data != {}:
            data = self.dict_data
        elif self.dict_data != data:
            self.dict_data = data
        print('\nFinding the optimal number of subspaces\nPROGRESS:')
        for subspace in range(subspace_min, subspace_max+1): # LOOP OVER THE RANGE OF SUBSPACES THAT IS GIVEN
            self.subspaces = subspace
            self.clustering()
            
            if measure == 'Sil':
                try: # CALCULATE THE SILHOUETTE SCORE AND SAVE IT
                    self.silhouette_score()
                    self.dict_silscores[subspace] = self.Sil_score
                except: # WHEN GIVEN NaN (WHEN NOT ALL CLUSTERS ARE FILLED OR ONLY 1 CLUSTER EXISTS) RETURN A VALUE THAT IS OUT OF THE SILHOUETTE RANGE
                    self.dict_silscores[subspace] = -2
  
            else: # If not 'Sil', the measure is 'E' (see assert commands above)
                self.Escore()
                self.dict_Escores[subspace] = self.E_score
                    
            print(f'{(((subspace-subspace_min)/(subspace_max-subspace_min))*100):.2f}%') # PRINT AT EVERY ITERATION HOW FAR THE EVALUATION PROCES IS

        if measure == 'Sil':
            self.subspaces =  list(self.dict_silscores.keys())[np.argmax(list(self.dict_silscores.values()))]
            self.Sil_score = self.dict_silscores[self.subspaces]
            self.clustering()
            self.Escore()
        else:
            self.subspaces = list(self.dict_Escores.keys())[np.argmin(list(self.dict_Escores.values()))]
            self.E_score = self.dict_Escores[self.subspaces]
            self.clustering()
            self.silhouette_score()

        return self.dict_clustered


"""
    PRINT TEXT FILES FORMATTED/SORTED LIKE THE INPUT FILE
    
"""
def return_txt_file(data, name, format_data):
    dict_unfolded = dict()
    score = False
    anticlusters = []
    for cluster in data:
        if cluster == -1:
            for i in data[cluster]:
                anticlusters.append(data[cluster][i])
        if type(data[cluster]) == int: 
            dict_unfolded[data[cluster]] = cluster
        else:
            try:
                for ID in data[cluster]:
                    dict_unfolded[ID] = cluster
            except:
                score = True
                break
    
    file = open(name, 'w')
    
    if score == False:
        for INDEX in format_data:
            try:
                line = str(INDEX) + " " + str(dict_unfolded[INDEX]) + "\n"
                file.write(line)
            except Exception:
                pass
    else:
        for subspace in data:
            line = str(subspace) + " " + str(data[subspace]) + "\n"            
            file.write(line)
            
    file.close()
    
    if len(anticlusters) > 0:
        print(f'\nFunction: {name} --> These clusters are lost (size is below threshold): {anticlusters}')
    

def main(filename):
     # APPLY KMCA TO THE GIVEN DATASET:
     dict_data = file_to_dict(filename)
     algorithm = str(input('Which clustering algorithm do you want to use?\n"KMCA" or "GBCA":\n '))
     if algorithm == 'KMCA':
         k_min = int(input('Minimum number of clusters (k_min):\n '))
         k_max = int(input('Maximum number of clusters (k_max):\n '))
         measure = str(input('By which score should the quality of clustering be assessed?\nType "E" or "Sil":\n '))
         seeds = list(range(int(input('How many different seeds should be chosen? (used for random initial cluster assignment)\nNumber of different seeds:\n '))))
         
         kmca = KMCA(data=dict_data)
         kmca_results = kmca.optimize(k_min=k_min, k_max=k_max, measure=measure, seeds=seeds)
         outfile_name = 'kmca_results.txt'
         return_txt_file(kmca_results, outfile_name, dict_data)
         print(f'\nClustering finished and optimized; The optimal number of clusters is {kmca.k}\nThe results of the best clustering are written to {outfile_name}\nThe E-score is {kmca.E_score}\nThe Silhouette score is {kmca.Sil_score}\n')
         
     elif algorithm == 'GBCA':
         scope = int(input('What is the range (in subspaces) in which a cell must be to another cell to be considered its neighbour?\n '))
         thres = int(input('What is the minimum amount of datapoints that must be in a cluster to consider it a valid cluster?\n '))
         subspace_min = int(input('Minimum number of subspaces (subspace_min):\n '))
         subspace_max = int(input('Maximum number of subspaces (subspace_max):\n '))
         measure = str(input('By which score should the quality of clustering be assessed?\nType "E" or "Sil":\n '))

         gbca = GBCA(data=dict_data, scope=scope, thres=thres)
         gbca_results = gbca.optimize(subspace_min=subspace_min, subspace_max=subspace_max, measure=measure)
         outfile_name = 'gbca_results.txt'
         return_txt_file(gbca_results, outfile_name, dict_data)
         print(f'\nClustering finished and optimized; The optimal number of subspaces is {gbca.subspaces}\nThe results of the best clustering are written to {outfile_name}\nThe E-score is {gbca.E_score}\nThe Silhouette score is {gbca.Sil_score}\n')
     else:
         raise ValueError('Value Error: no valid algorithm name provided. Please type either "KMCA" or "GBCA".') 
     return outfile_name