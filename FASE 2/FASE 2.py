import numpy as np

"""

   DEFINITION TO IMPORT TEXTFILES AND RETURN DICTIONARIES:
        
"""

def file_to_lib(FileName):
    def openFile(FileName):
        read_data = open(f'{FileName}')
        data = read_data.read()
        data = data.splitlines()
        read_data.close()
        
        return data
        
    lib = {}
    
    for lines in openFile(FileName):
        line = lines.split()
        if len(line[1:]) == 1:
            lib[int(line[0])] = int(line[1])
        else:
            lib[int(line[0])] = line[1:]
            
    return lib



"""

    K-MEANS CLUSTERING ALGORITHM (KMCA):
        
"""

class KMCA:
    def __init__(self, k=6, seed=20):
        '''
        preconditions:  The number of clusters for, the number of seed for seed, the maximum number of times
                        you want to itterate over the data
        postconditions: sets the number of clusters to self.k
                        sets the number of seeds to self.seed
        '''
        self.k = k + 1
        self.seed = seed
    
    def fit(self, data):
        '''
        preconditions:  a library with that contains the data. the library must have tehe gene code as key and a list
                        of the datapoints as value.
        postconditions: returns a library: 
                            keys   --> clusters
                            values --> genes that are classefied in that cluster
        this functions task: to assign all genes to its right cluster.
        '''
        
        self.input_data = data
        
        def E_score():
            '''
            preconditions:  /
            postconditions: returns the e_score as an integer.
            this functions task: calculate the e_score for the clusterd genes 
            '''
            
            mean = 0

            for k in range(1,self.k):
                score = np.sum([np.square(abs(np.subtract(self.centroids[k], self.lib_norm[dis]))) for dis in self.labeled_data[k]])
                mean += score

            e_score = mean / self.k

            return e_score

        def normalize():
            '''
            preconditions:  /
            postconditions: sets the normalised coordinates as content for the variable self.lib_norm
            this functions task: normalise the coordinates of all genes
            '''
            
            normalized = dict()
            for i in data:
                norm = np.sqrt(np.sum(np.square(np.array(data[i]).astype(float)))) ## SET DATA TYPE TO FLOAT --> THEN CALCULATE THE LENGTH OF THE VECTOR BY APPLYING THE FORMULA IN THE CASUS
                normalized[i] = np.divide(np.array(data[i]).astype(float), norm) ## DIVIDE EACH VALUE IN THE CORDINATES BY THE LENGTH OF THE VECTOR
            
            self.lib_norm = normalized
        
        
        def rand_label():
            '''
            preconditions:  /
            postconditions: returs a library:
                                keys   --> cluster
                                values --> genes that are random classefied to that cluster
            this functions task: random assigning genes to a cluster
            '''
            
            labeld_data = dict()
            labels_id = dict()
            np.random.seed(self.seed) # SET THE RANDOMNESS OF NUMPY
            
            # MAKE A DUMMY LIST WITH THE LENGTH OF THE INPUT DATA WITH EACH A RANDOM VALUE IN THE RANGE OF K
            labels = np.random.randint(1,self.k,len(self.lib_norm)) 

            # MAKE DICTIONARY WITH AS KEY THE PROTEIN-ID AND AS VALUE THE RANDOM CLUSTER-VALUE FROM THE DUMMY
            for i in range(len(self.lib_norm)):
                labeld_data[list(self.lib_norm.keys())[i]] = labels[i] 

            # MAKE DICTIONARY WHERE THE KEY IS THE CLUSTER-#NR AND THE VALUE IS A LIST OF PROTEIN-IDs THAT ARE IN THAT CLUSTER
            for K in range(1,self.k):
                labels_id[K] = np.array([i for i,j in labeld_data.items() if j == K]) 
            
            return labels_id
            
        def centroid():
            '''
            preconditions:  /
            postconditions: returns a library:
                                keys   --> cluster
                                values --> list of coordinats for the center of that cluster
            this functions task: calculating the centers of the clusters
            '''
            
            centroids = dict()
            
            # CALCULATE THE MEAN VECTOR PER CLUSTER AND SAFE IT AS CENTROID OF THE CLUSTER
            for k in range(1,self.k):
                centroids[k] = np.mean(np.array([self.lib_norm[index] for index in self.labeled_data[k]]), axis=0) 

            return centroids
        
        def assign_cluster(unlabeld_data):
            '''
            preconditions:  a library:
                                keys  --> gen
                                value --> the normalised coordinates of the gen
            postconditions: a library:
                                keys  --> cluster
                                value --> genes classified to that cluster                     
            this functions task: assigning genes to the cluster with the nearest center
            '''
            
            lib_clust = {k:[] for k in range(1,self.k)}
            
            for cords in unlabeld_data:
                # GET THE CLUSTER NUMBER WHERE THE DISTANCE IS THE SMALLEST FROM VECTOR TO CENTROID
                # np.armin --> GIVES THE INDEX OF THE LOWEST VALUE IN A LIST. THIS RANGES FROM 0 to k-1, SO 1 NEEDS TO BE ADDED AT THE END
                clust = np.argmin([np.sqrt(np.sum(np.square(np.subtract(self.centroids[k],unlabeld_data[cords])))) for k in range(1,self.k)])+1
                lib_clust[clust].append(cords)

            return lib_clust
        
        # CALL ALL FUNCTIONS FOR THE FIRST TIME SO THE DATA HAS A STARTING POINT
        normalize()
        self.labeled_data = rand_label()
        self.centroids = centroid()
        self.E_score = E_score()
        
        # CHECKING VAR TO LOOP THE ALGORITHM UNTIL THE BEST CENTRIODS ARE MADE
        optimized = False
        
        while optimized == False:
            self.labeled_data = assign_cluster(self.lib_norm)

            self.centroids = centroid()
            new_e_score = E_score()

            # CHECK IF THERE IS AN IMPROVEMENT IN THE E-SCORE --> IF NOT THEN THE BEST CENTROIDS ARE PICKED AND THE OPTIMISATION LOOP NEEDS TO BREAK
            if self.E_score > new_e_score:
                self.E_score = new_e_score
            else:
                optimized = True
                    
        return self.labeled_data
    
    def optimize(self, min_seed=0, max_seeds=30):
        '''
        preconditions:  min,max seed: range to loop the seeds in
        postconditions: a library:
                            keys  --> cluster
                            value --> genes classified to that cluster                     
        this functions task: assigning genes to the cluster with the best centroids
        '''
        
        # START WITH THE BIGGEST POSSIBLE #NR SO THE FIRST TESTED SEED IS ALWAYS LOWER THAN THIS #NR SO THE IF STATEMENT WORKS CORRECTLY
        old_e = float('inf') 
        
        # TEST THE SEEDS FOR WHICH THE RANGE IS SET:
        for seed in range(min_seed,max_seeds+1):
            self.seed = seed
            self.fit(self.input_data)
            new_e = self.E_score
            
            # COMPARE THE NEWLY GENERATED E-SCORE WITH THE OLD-BEST E-SCORE --> IF NEWLY GENERATED E-SCORE IS BETTER THEN THE LAST ONE SAFE IT
            if new_e < old_e: 
                best_labeld_data = self.labeled_data
                best_seed = seed
                old_e = new_e
                
        self.seed = best_seed
        
        return best_labeld_data, best_seed, old_e
        
    
    def predict(self, input_data):
        lib_clust = {k:[] for k in range(1,self.k)}
        
        # CHECK IF INPUT TYPE CONTAINS MULTIPLE DATA POINTS OR JUST A SINGLE DATA POINT:
        if type(input_data) == dict:
            for cords in input_data:
                # NORMALIZE INPUT DATA
                norm = np.sqrt(np.sum(np.square(np.array(input_data[cords]).astype(float)))) 
                normalized = np.divide(np.array(input_data[cords]).astype(float), norm)

                # CALCULATE CLOSEST CENTROID AND SAVE IT TO CORRESPONDING CLUSTER
                clust = np.argmin([np.sqrt(np.sum(np.square(np.subtract(self.centroids[k],normalized)))) for k in range(1,self.k)])+1
                lib_clust[clust].append(cords)
                
        if type(input_data) == list:
            # NORMALIZE INPUT DATA
            norm = np.sqrt(np.sum(np.square(np.array(input_data).astype(float))))
            normalized = np.divide(np.array(input_data).astype(float), norm)
            
            # CALCULATE CLOSEST CENTROID AND SAVE IT TO CORRESPONDING CLUSTER
            clust = np.argmin([np.sqrt(np.sum(np.square(np.subtract(self.centroids[k],normalized)))) for k in range(1,self.k)])+1
            lib_clust[clust].append(normalized)
            
        return lib_clust
    


"""

    OUR OWN CLUSTER METHOD: GRID based CLUSTERING ALGORITHM (GCA)
    
"""

class GCA:    
    def __init__(self, r=1, j=10, tresh=4):
        '''
        preconditions:  The range of neigbors (r), the total amount of sub-spaces per axis (j), the treshhold on how many 
                        points there need to be per cel to count it as clustercel (tresh).
        postconditions: sets the range to self.r
                        sets the amount of sub-spaces to self.j
                        sets the treshhold to self.tresh
        '''
        self.r = r
        self.j = j
        self.tresh = tresh    
    
    def fit(self, data):
        '''
        preconditions:  a library with that contains the data. the library must have the gene code as key and a list
                        of the datapoints as value.
        postconditions: returns a library: 
                            keys   --> clusters
                            values --> unique cell ID
        this functions task: to assign all genes to a cell.
        '''
        self.indexes = list(data.keys())

        def axis_data(data):
            '''
            preconditions:  a library with that contains the data. the library must have tehe gene code as key and a list
                            of the datapoints as value.
            postconditions: returns a library: 
                                keys   --> axis
                                values --> list of all values in that specific axis
            this functions task: returning a library with as key every axis and as value a tuple of all value of the axis.
            '''
            cordinates = dict()

            for cordinate in range(len(data[self.indexes[0]])): # GET THE GENERAL LENGTH OF THE CORDINATES --> SO FOR EVERY AXIS THERE WILL BE A LIST
                x = list()
                for index in data: # GET EVERY PROTEIN ID IN THE DATA DICTIONARY
                    x.append(float(data[index][cordinate]))  # APPEND THE AXIS LIST WITH THE CORROSPONDING AXIS FROM EVERY PROTEIN

                cordinates[cordinate] = np.array(x) # MAKE TUPLE OF THE AXIS VALUES
            return cordinates

        def interval_data(axis_data):
            '''
            preconditions:  a library that contains the axis-data.
            postconditions: return per axis the intervals of the subspaces
            this functions task: returning a library with the subspaces of the data sorted on the axis.
            '''
            dic = dict()

            # GET THE RANGE (MIN-MAX) OF EVERY AXIS --> SO SUB-CELLS CAN BE MADE IN THIS RANGE
            for axis in axis_data:
                minim = np.min(axis_data[axis])
                maxim = np.max(axis_data[axis])

                dic[axis] = np.mgrid[minim:maxim:complex(0,self.j)] # J IS THE #NR OF INTERVALLS WE WANT EVERY AXIS

            return dic

        def apply_grid(data):
            '''
            preconditions:  - a library with that contains the data. the library must have the gene code as key and a list
                              of the datapoints as value. 
                            
            postconditions: returns a library: 
                                keys   --> gene ID
                                values --> unique cel ID
            this functions task: returning a library with gene ID as key and as value the unique cel ID.
            '''
            
            grid = list()
            self.lib_cells = dict()
            
            lib_cords = axis_data(data) # GET A DICTIONARY WITH ALL THE DATA PER AXIS
            lib_interv = interval_data(lib_cords) # GET A DICTIONARY WITH ALL THE INTERVALS PER AXIS
            
            # GET FOR EVERY AXIS THE SUBCEL-ID
            for axis in lib_cords: # LOOP OVER EVERY AXIS
                for subspace in range(1,len(lib_interv[axis])): # LOOP OVER EVERY SUBSPACE (SPACE BETWEEN 2 INTERVALS)
                    # CONDITIONS ARE APPLIED TO THE VALUES OF AN AXIS
                    # --> CONDITIONS ARE THAT THE VALUES NEED TO BE BETWEEN 2 INTERVALS SET BY THE INTERVAL_DATA FUNCTION
                    # --> LIB_CORDS => A LIST WITH AS VALUES THE SUBSPACE-ID 
                    lib_cords[axis][(lib_cords[axis] > lib_interv[axis][subspace-1]) & (lib_cords[axis] < lib_interv[axis][subspace])] = subspace
                
                # APPEND THE GRID LIST WITH THE CALCULATED SUBSPACE LISTS PER AXIS
                grid.append(lib_cords[axis])
            
            # CONVERT THE GRID LIST TO A DICTIONARY
            for i in range(len(self.indexes)):
                self.lib_cells[self.indexes[i]] = list()
                
                for axis in grid:
                    self.lib_cells[self.indexes[i]].append(int(axis[i]))
                    
            return self.lib_cells
                
        return apply_grid(data)
        
    def clustering(self, cell_ids):
        '''
        preconditions:  a library with that contains the data. the library must have the gene code as key and a list
                        of the datapoints as value.
        postconditions: returns a library: 
                            keys   --> clusters
                            values --> genes that are classefied in that cluster
        this functions task: to assign all genes to its right cluster.
        '''
        def get_neigbours(ID):
            '''
            preconditions:  a cell ID.  
            postconditions: returns a list with all the IDs of the neigboring cells
            this functions task: returning neigbor cel IDs.
            '''
            neigb = list()
            #GET EVERY NEIGBOR OF THE INPUT-CELL:
            for indx in cell_ids: # LOOP OVER EVERY CELL-INDEX THAT IS GENERATED
                # THE DISTANCE TO THE INPUT-CELL IS A FORM OF MANHATTAN DISTANCE --> SUBTRACT EVERY SUBSPACE-ID FROM THE INPUT-CELL SUBSPACES
                # --> WHEN THE OUTCOME OF EVERY SUBSPACE OF A CELL-ID IS IN THE SET RANGE, IT CAN BE SEEN AS A NEIGBOR
                distance = np.subtract(self.lib_cells[ID],self.lib_cells[indx],out=np.array([0 for i in range(8)])) 
                if all([-self.r < i < self.r for i in distance]): # CHECK IF ALL THE VALUES OF THE SUBSTRACTION (in distance) ARE IN RANGE
                    neigb.append(indx) # APPEND THE NEIGBOR LIST WITH THE PROTEIN ID OF THE NEIGBOR

            neigb.remove(ID)
            
            return neigb
    
        lib_cluster = dict()
        
        def get_all_neigbours(ID,location):
            '''
            preconditions:  a cell ID + a location in the list of cell-neigbors
            postconditions: updated location + updated group (all neigbors) list
            this functions task: gathering all neigbors of a originated cell
            '''
            
            self.group += get_neigbours(ID)
            self.group = list(np.unique(self.group))
            
            # UPDATE LOCATION
            location += 1
            
            # IF THE LOCATION IS GREATER THAN THE LIST OF ALL NEIGBOURS --> THE LIST OF NEIGBOURS HAS ALL ITS NEIGBOURS ASSIGNED TO IT
            # IF THIS IS NOT YET THE CASE CALL THE SAME FUNCTION WITH THE UPDATED LOCATION
            if len(self.group) > location:
                get_all_neigbours(self.group[location], location)
        
        def make_clusters(indexes, clust_nr=1):
            '''
            preconditions:  list of all indexes 
            postconditions: updated list of indexes that are not been assigned to a neigbor group
            this functions task: - calling the function get_all_neigbours efficiently so no duplicates are formed --> get_all_neigbours wont form duplicate neigbour groups
                                 - returning a dictionary with every cluster as key and the corrisponding protein-IDs in a list as value  
            '''
            if len(indexes) > 0:
                indexes = indexes.copy()
                
                self.group = list() # FOR EVERY PROTEIN INDEX MAKE A NEW GROUP
                get_all_neigbours(indexes[0], -1) # CALL THE FUNCTION GET ALL NEIGBOURS FOR THE FIRST INDEX IN THE LIST OF INDEXES, STARTING POINT == -1, BECAUSE -1 IS UPDATED BY 1 => STARTING POINT IS REALY 0
                
                if len(self.group) == 0:
                    # IF ORIGINAL PROEIN INDEX HAS NO NEIGBOURS, THEN GIVE THE GROUP ONLY THE ONE PROTEIN-ID AND REMOVE THAT ONE ID FROM THE LIST OF INDEXES
                    cluster_group = indexes[0]
                    indexes.remove(cluster_group)
                else:    
                    # ELSE: GET EVERY INDEX IN THE GROUP LIST AND REMOVE IT FROM THE LIST OF INDEXES
                    cluster_group = self.group
                    for index in cluster_group:
                        try:
                            indexes.remove(int(index))
                        except:
                            pass
                
                # MAKE CLUSTER DICTIONARY
                lib_cluster[clust_nr] = cluster_group
                # UPDATE CLUSTER #NR BY 1 EVERY ITERATION
                clust_nr += 1
                # CALL THE SAME FUNCTION WITH THE UPDATED LIST OF INDEXES + THE UPDATED CLUSTER #NR
                make_clusters(indexes, clust_nr)
                
            return lib_cluster
        self.lib_cluster = make_clusters(self.indexes, 1)
        
        return self.lib_cluster
"""

    GET RESULTS BY CALLING CLASSES
    
"""
        
lib_data = file_to_lib('Data\Voorbeeld_clusterdata.txt')
lib_results = file_to_lib('Data\Voorbeeld_clusterresult.txt')

kmca = KMCA()
kmca.fit(lib_data)
kmca_result, kmca_best_seed, kmca_e_score = kmca.optimize(10,30)

gca = GCA(r=1, j=5)
gca_fit_data = gca.fit(lib_data)
gca_results = gca.clustering(gca_fit_data)

"""

    PRINT TEXT FILES FORMATTED/SORTED LIKE THE INPUT FILE
    
"""

def return_txt_file(data, name, format_data=lib_data):
    lib_unfolded = dict()
    
    for cluster in data:
        if type(data[cluster]) == int:
            lib_unfolded[data[cluster]] = cluster
        else:
            for ID in data[cluster]:
                lib_unfolded[ID] = cluster
    
    file = open(f'{name}results.txt', 'w')
    indexes_lost = list()
    
    for INDEX in format_data:
        try:
            line = str(INDEX) + " " + str(lib_unfolded[INDEX]) + "\n"
            file.write(line)
        except:
            indexes_lost.append(INDEX)
    
    print(f'Function: {name} --> these indexes are lost: {indexes_lost}')
        
    file.close()
    
return_txt_file(kmca_result, 'kmca_')
return_txt_file(gca_results, 'gca_')