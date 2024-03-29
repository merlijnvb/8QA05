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
   GRID BASED CLUSTERING:
        
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
                        --> self.lib_[...] = empty dictionary that is filled later.
                            --> lib_data can be an empty dicionary or can be filled from the beginning. This choice is up to the user.
                        --> self.indices_not_clustered = empty array that is filled later.
                        --> self.E_score = set to infinity, because if we don't cluster the data in the first iteration the E_score is not representative.
                        --> self.Sil_score = set to 0, because if we don't cluster the data in the first iteration the Sil_score is not representative.
                        --> self.lib_silscores = empty dictionary that is later filled with different Silhouette scores when parameters are changed.
                        --> self.lib_Escores = empty dictionary that is later filled with different E-scores when parameters are changed.
        '''
        self.scope = scope
        self.subspaces = subspaces 
        self.thres = thres 
        self.lib_data = data
        self.lib_axis = {}
        self.lib_interval = {}
        self.lib_grid = {}
        self.lib_clustered = {}    
        self.indices_not_clustered = np.array([])
        self.E_score = float('inf')
        self.Sil_score = -2
        self.lib_Escores = {}
        self.lib_silscores = {}
            
        
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
        assert data != {}, "No input given or empty dictionary given as input for 'data'. --> interval_data(data)"
        assert type(data) == dict, "Input of 'data' is of a type other than dictionary, must be dictionary. --> interval_data(data)"
        if data != self.lib_data:
            self.lib_data = data
        
        # CONVERT THE DATA IN LIB_DATA TO A DICTIONARY WITH THE KEYS BEING THE AXES AND THE VALUES BEING THE CORRESPONDING VALUES OF EVERY AXIS FOR EVERY PROTEIN
        for axis in range(len(list(self.lib_data.values())[0])): 
            x = []
            for index in self.lib_data: 
                x.append(float(self.lib_data[index][axis]))  
        
            self.lib_axis[axis] = np.array(x) 
       
        # CALCULATE THE INTERVALS PER AXIS USING THE MINIMUM AND MAXIMUM (RANGE IN WHICH ALL DATA POINTS LAY)
        for axis in self.lib_axis:
            minimum = np.min(self.lib_axis[axis])
            maximum = np.max(self.lib_axis[axis])

            self.lib_interval[axis] = np.mgrid[minimum:maximum:complex(0,self.subspaces+1)] # MAKE 1-D GRID => ((abs(min)+max/subspaces) * n_step) + min
        
        
        
    def grid_data(self):
        '''
        Postconditions: A dictionary:
                        --> Key = protein ID
                        --> Value = list of integers which represent the subspaces of the axes (the list represents the cell ID)
        Task of function: calculate the unique cell ID for every protein
        '''
        
        grid = {}
        
        for axis in self.lib_interval: # LOOP OVER EVERY AXIS
            lib_axis_subspace = self.lib_axis[axis].copy() # COPY THE AXIS LIST TO CONVERT THE DATAPOINTS WHICH FALLS BETWEEN INTERVALS TO SUBSPACE IDs
            
            for interval in range(1,len(self.lib_interval[axis])): # LOOP OVER EVERY INTERVAL EXCEPT THE FIRST ONE
                # APPLY THE INTERVAL CONDITIONS TO THE CORRESPONDING AXIS AND CONVERT THE DATAPOINTS WHICH RETURN TRUE TO THE CONDITIONS TO THE SUBSPACE ID
                lib_axis_subspace[(self.lib_axis[axis] >= self.lib_interval[axis][interval-1]) & (self.lib_axis[axis] <= self.lib_interval[axis][interval])] = interval
                         
            # SAVE A LIBRARY WITH AS KEY THE AXIS AND AS VALUE A LIST OF THE ASSIGNED SUBSPACE IDs FOR EVERY PROTEIN
            grid[axis] = lib_axis_subspace.astype(int)
            
        for i in range(len(self.lib_data)): # LOOP OVER EVERY PROTEIN
            # CONVERT THE GRID LIBRARY TO A LIBRARY WITH THE KEYS SET TO THE PROTEIN ID AND THE VALUE SET TO A LIST CONTAINING ALL THE SUBSPACE IDs FOR EVERY AXIS ==> CELL ID
            index = list(self.lib_data.keys())[i]
            self.lib_grid[index] = np.zeros((len(self.lib_interval),),dtype=int)
            
            for axis in grid:
                self.lib_grid[index][axis] = grid[axis][i]        
    
    
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
            if all([-self.scope <= axis_distance <= self.scope for axis_distance in np.subtract(self.lib_grid[index], self.lib_grid[index0])]):
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
        assert data != {}, "No input given or empty dictionary given as input for 'data'. --> optimize(data) or clustering(data) or __init__(data) or interval_data(data)"
        assert type(data) == dict, "Input of 'data' is of a type other than dictionary, must be dictionary. --> clustering(data)"
        if self.lib_data != data:
            self.lib_data = data
        
        self.indices_not_clustered = np.array(list(self.lib_data.keys()))
        self.interval_data(data)
        self.grid_data()
        
        cluster_nr = 0
        anticluster_nr = 0
        lib_anticlustered = {}
        
        for index in self.indices_not_clustered: # LOOP OVER EVERY PROTEIN THAT IS NOT YET ASSIGNED TO ANY CLUSTER OR ANTICLUSTER
            if index in self.indices_not_clustered: 
                cluster = []
                self.assign_neigbours(index,cluster)
                
                # CHECK IF THE CLUSTER IS BIGGER THEN THE SET THRESHOLD
                # --> IF TRUE: ADD THE CLUSTER TO THE LIBRARY CONTAINING THE CLUSTERS
                # --> IF FALSE: ADD THE CLUSTER TO THE LIBRARY CONTAINING THE NON-CLUSTERS IN THE LIBRARY CONTAINING THE CLUSTERS
                if len(cluster) > self.thres:              
                    self.lib_clustered[cluster_nr] = cluster
                    cluster_nr += 1
                else:
                    lib_anticlustered[anticluster_nr] = cluster
                    anticluster_nr += 1
            
        # APPEND THE CLUSTER LIBRARY WITH THE LIBRARY CONTAINING THE NON-CLUSTERS             
        self.lib_clustered[-1] = lib_anticlustered
        
        return self.lib_clustered
    
    
    
    def Escore(self):
        '''
        Postconditions: E_score => score how good the grid-based clustering fit is
        Task of function: calculation the goodness of the k-means clustering application of certain centroids
        '''
        
        summation = 0
        lib_centroid = {}
        
        # CALCULATE THE CENTROIDS OF EVERY CLUSTER
        for k in self.lib_clustered: # LOOP OVER THE CLUSTERS
            if k >= 0: # ONLY CALCULATE THE CENTROIDS FOR THE CLUSTERS NOT FOR THE NOT ASSIGNED DATAPOINTS
                coordinate = []
                for axis in self.lib_axis: # LOOP OVER THE AXES
                    mean_x = 0
                    
                    for protein in [value for index, value in self.lib_data.items() if index in self.lib_clustered[k]]: # CALCULATE THE MEAN VALUE OF THE AXIS
                        mean_x += protein[axis]
                    
                    coordinate.append(mean_x / len(self.lib_clustered[k]))
                    lib_centroid[k] = coordinate
        
        for k in lib_centroid: # CALCULATE FOR EVERY CLUSTER THE E_SCORE
            summation += np.sum([np.square(abs(np.subtract(lib_centroid[k], self.lib_data[dis]))) for dis in self.lib_clustered[k]])
            
        # CALCULATE THE MEAN E_SCORE FOR THE ENTIRE GRID CLUSTERING FIT
        self.E_score = summation / len(lib_centroid)
        
        return
    
    
    
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
        
        for index in self.lib_data:
            # GET THE CLUSTER #NR WHICH THE PROTEIN IS ASSIGNED TO
            cluster_nr = [cluster for cluster, indices in self.lib_clustered.items() if (index in indices) and (cluster != -1)] # GET CLUSTER #NR THE PROTEIN IS ASSIGNED TO
            
            #CHECK IF THE PROTEIN IS ASSIGNED TO A CLUSTER
            if cluster_nr != []:
                cluster_nr = cluster_nr[0]
            
                if len(self.lib_clustered[cluster_nr]) > 1:
                    # CALCULATE THE MEAN DISTANCE OF THE PROTEIN TO ALL THE OTHER PROTEINS IN ITS CLUSTER (INTERNAL):
                    dissimilarity_internal = np.sum([np.linalg.norm(self.lib_data[index] - self.lib_data[index_in]) for index_in in self.lib_data]) / (len(self.lib_data)-1)
                    
                    # CALCULATE THE MINIMAL MEAN DISTANCE OF THE PROTEIN TO ALL THE OTHER PROTEINS FROM THE OTHER CLUSTERS (EXTERNAL):
                    # --> MEAN DISTANCE IS GROUPED PER CLUSTER. FROM THIS LIST THE MINIMAL DISTANCE IS CALCULATED
                    dissimilarity_external = np.min([np.mean([np.linalg.norm(self.lib_data[index] - self.lib_data[index_ex]) for index_ex in self.lib_data if index_ex in self.lib_clustered[cluster]]) for cluster in self.lib_clustered if (cluster != cluster_nr) & (cluster != -1)])
                    
                    # CALCULATING THE SILHOUETTE SCORE FOR PROTEIN BY APPLYING THE FORMULA:
                    sil_score_i = (dissimilarity_external - dissimilarity_internal) / max(dissimilarity_external,dissimilarity_internal)
                    
                if len(self.lib_clustered[cluster_nr]) == 1:
                    sil_score_i = 0
                
                list_scores.append(sil_score_i)
            
        self.Sil_score = np.mean(list_scores)
        
        if np.isnan(self.Sil_score):
            self.Sil_score = float('-inf')
                
        return
    
    def optimize(self, data={}, subspace_min=10, subspace_max=20, measure = 'Sil'):
        '''
        Preconditions:  data: a dictionary:
                            --> Key = protein ID
                            --> Value = list of floats
                        subspace_min, subspace_max: range in which parameter subspaces can be varied
                        measure: measure by which to evaluate how well the clustering fits (either 'E' for E-score or 'Sil' for silhouette score)
        Postconditions: lib_clustered: a dictionary:
                                    --> Key = Cluster #nr
                                    --> Value = list of protein IDs that are assigned to that cluster when the E_score (the goodness of the fit) is in its maximum
                        lib_Silscores/lib_Escores: a dictionary:
                                                --> Key = subspace
                                                --> Value = Sil_score/E_score
        Task of function: calculating which subspace gives the best results --> returning the clustered proteins for which the Sil_score is at its highest when comparing different subspaces
        '''
        assert data != {}, "No input given or empty dictionary given as input for 'data'. --> optimize(data) or clustering(data) or __init__(data) or interval_data(data)"
        assert type(data) == dict, "Input of 'data' is of a type other than dictionary, must be dictionary. --> clustering(data)"
        if self.lib_data != data:
            self.lib_data = data
        assert subspace_min < subspace_max, "subspace_min must be smaller than subspace_max. --> optimize(subspace_min,subspace_max)"
        assert (measure == 'Sil' or measure == 'E'), "Measure must either be 'Sil' or 'E'. --> optimize(measure)"
        
        for subspace in range(subspace_min, subspace_max+1): # LOOP OVER THE RANGE OF SUBSPACES THAT IS GIVEN
            self.subspaces = subspace
            self.clustering(data)
            
            if measure == 'Sil':
                try: # CALCULATE THE SILHOUETTE SCORE AND SAVE IT
                    self.silhouette_score()
                    self.lib_silscores[subspace] = self.Sil_score
                except: # WHEN GIVEN NaN (WHEN NOT ALL CLUSTERS ARE FILLED OR ONLY 1 CLUSTER EXISTS) RETURN A VALUE THAT IS OUT OF THE SILHOUETTE RANGE
                    self.lib_silscores[subspace] = -2
                    
            else: # If not 'Sil', the measure is 'E' (see assert commands above)
                self.Escore()
                self.lib_Escores[subspace] = self.E_score
                    
            print(f'{((subspace-subspace_min)/(subspace_max-subspace_min))*100}% ', end='\r') # PRINT AT EVERY ITERATION HOW FAR THE EVALUATION PROCES IS
        
        if measure == 'Sil':
            self.subspaces =  list(self.lib_silscores.keys())[np.argmax(list(self.lib_silscores.values()))]
        else:
            self.subspaces = list(self.lib_Escores.keys())[np.argmax(list(self.lib_Escores.values()))]
            
        self.clustering(data)
        self.Escore()
        self.silhouette_score()
        
        return self.lib_clustered


lib_data = file_to_lib('Data\Voorbeeld_clusterdata.txt')
lib_results = file_to_lib('Data\Voorbeeld_clusterresult.txt')
   
gbca = GBCA(data=lib_data)
gbca_results = gbca.optimize(lib_data,8,9,'Sil')
gbca_scores = gbca.lib_Escores


def return_txt_file(data, name, format_data=lib_data):
    lib_unfolded = dict()
    score = False
    anticlusters = []
    for cluster in data:
        if cluster == -1:
            for i in data[cluster]:
                anticlusters.append(data[cluster][i])
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
    
    if score == False:
        for INDEX in format_data:
            try:
                line = str(INDEX) + " " + str(lib_unfolded[INDEX]) + "\n"
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
    
return_txt_file(gbca_results, 'gca_')
return_txt_file(gbca_scores, 'gca_e_')
