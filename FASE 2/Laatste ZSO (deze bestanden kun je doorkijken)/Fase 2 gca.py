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
    def __init__(self, scope=1, subspaces=9, thres=3, data={}):
        self.scope = scope
        self.subspaces = subspaces 
        self.thres = thres 
        self.lib_data = data
        self.lib_axis = {}
        self.lib_interval = {}
        self.lib_cells = {}
        self.lib_grid = {}
        self.lib_clustered = {}
        self.indexes_not_clusterd = np.array(list(self.lib_data.keys()))
        self.E_score = 0
        self.Sil_score = 0
        
    def interval_data(self, data={}):
        if self.lib_data == {}:
            self.lib_data = data
        
        
        
        for axis in range(len(list(self.lib_data.values())[0])): 
            x = []
            for index in self.lib_data: 
                x.append(float(self.lib_data[index][axis]))  
        
            self.lib_axis[axis] = np.array(x) 
       
        for axis in self.lib_axis:
            minimum = np.min(self.lib_axis[axis])
            maximum = np.max(self.lib_axis[axis])

            self.lib_interval[axis] = np.mgrid[minimum:maximum:complex(0,self.subspaces+1)]
        
    def grid_data(self):
        grid = {}
        
        for axis in self.lib_interval:
            lib_axis_subspace = self.lib_axis[axis].copy()
            
            for interval in range(1,len(self.lib_interval[axis])):                
                lib_axis_subspace[(self.lib_axis[axis] >= self.lib_interval[axis][interval-1]) & (self.lib_axis[axis] <= self.lib_interval[axis][interval])] = interval
                         
                
            grid[axis] = lib_axis_subspace.astype(int)
            
        for i in range(len(self.lib_data)):
            index = list(self.lib_data.keys())[i]
            self.lib_grid[index] = np.zeros((len(self.lib_interval),),dtype=int)
            
            for axis in grid:
                self.lib_grid[index][axis] = grid[axis][i]        
        
    
    
    def assign_neigbours(self, index0, cluster):
        neighbours = []
    
        for index in self.indexes_not_clusterd:
            (np.subtract(self.lib_grid[index], self.lib_grid[index0]))
            if all([-self.scope <= axis_distance <= self.scope for axis_distance in np.subtract(self.lib_grid[index], self.lib_grid[index0])]):
                neighbours.append(index)
                cluster.append(index)
                self.indexes_not_clusterd = np.delete(self.indexes_not_clusterd, np.where(self.indexes_not_clusterd == index))          
        
        if len(neighbours) > 1:
             for neighbour in neighbours:
                 self.assign_neigbours(neighbour,cluster)
        
        
    
    def clustering(self, data={}):
        if self.lib_data == {}:
            if data == {}:
                raise ImportError("No input given. --> clustering(data) or __init__(data) or interval_data(data)")
            else:
                self.lib_data = data
                self.indexes_not_clusterd = np.array(list(self.lib_data.keys()))
        if (self.indexes_not_clusterd.size <= 1) & (self.lib_data != {}):
            self.indexes_not_clusterd = np.array(list(self.lib_data.keys()))
         
        self.interval_data()
        self.grid_data()
        
        cluster_nr = 0
        anticluster_nr = 0
        lib_anticlustered = {}
        
        for index in self.indexes_not_clusterd:
            if index in self.indexes_not_clusterd:
                cluster = []
                self.assign_neigbours(index,cluster)
                
                if len(cluster) > self.thres:                 
                    self.lib_clustered[cluster_nr] = cluster
                    cluster_nr += 1
                else:
                    lib_anticlustered[anticluster_nr] = cluster
                    anticluster_nr += 1
                    
        self.lib_clustered[-1] = lib_anticlustered
        
        return self.lib_clustered
    
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
            cluster_nr = [cluster for cluster, indexes in self.lib_clustered.items() if (index in indexes) and (cluster != -1)] # GET CLUSTER #NR THE PROTEIN IS ASSIGNED TO
            
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
                
        return self.Sil_score
    
    def Escore(self):
        '''
        Postconditions: E_score => score how good the k-means clustering fit is
        Task of function: calculation the goodness of the k-means clustering application of certain centroids
        '''
        
        summation = 0
        lib_centroid = {}
        
        ## HIER VERDER WERKEN
        for k in self.lib_clustered:
            if k >= 0:
                cordinate = []
                for axis in self.lib_axis:
                    mean_x = 0
                    
                    for protein in [value for index, value in self.lib_data.items() if index in self.lib_clustered[k]]:
                        mean_x += protein[axis]
                    
                    cordinate.append(mean_x / len(self.lib_clustered[k]))
                    lib_centroid[k] = cordinate
        
        for k in lib_centroid: # CALCULATE FOR EVERY CLUSTER THE E_SCORE
            summation += np.sum([np.square(abs(np.subtract(lib_centroid[k], self.lib_data[dis]))) for dis in self.lib_clustered[k]])
            
        # CALCULATE THE MEAN E_SCORE FOR THE ENTIRE K-MEANS FIT
        self.E_score = summation / len(lib_centroid)
        
        return self.E_score
        
    def optimize(self, subspace_min=10, subspace_max=40):
        Sil_score_best = -2 # SET E_SCORE TO INFINITY SO THE 0TH ITERATION HAS THE LEAST FAVORABLE SCORE
        
        for subspace in range(subspace_min, subspace_max+1): # LOOP OVER THE RANGE OF SEEDS THAT IS GIVEN
            self.subspaces = subspace
            self.clustering()
            self.silhouette_score()
            
            if self.Sil_score > Sil_score_best: # CHECK IF THE GIVEN SEED HAS A BETTER FIT THAN THE LAST BEST ONE, IFSO SAVE IT
                Sil_score_best = self.Sil_score
                best_subsapce = subspace
                
        self.subspaces = best_subsapce 
        self.clustering()
        self.E_score()
        self.silhouette_score()
        
        return self.lib_clustered

lib_data = file_to_lib('Data\Voorbeeld_clusterdata.txt')
lib_results = file_to_lib('Data\Voorbeeld_clusterresult.txt')
   
gca = GCA(data=lib_data)
results = gca.optimize()
print(gca.Sil_score)