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
        preconditions:  k is the number of clusters to be made by the algorithm;
                        seed is the seed number used by the 'random' module
        postconditions: sets the number of clusters to self.k
                        sets the number of seeds to self.seed
        '''
        self.k = k + 1
        self.seed = seed
    
    def fit(self, data):
        '''
        Assigns all genes to their right clusters.
        preconditions:  data is a library that contains the data. the library must have the gene code as key and a list
                        of the datapoints as value.
        postconditions: returns a library: 
                            keys   --> clusters
                            values --> genes that are assigned to that cluster
        '''
        
        self.input_data = data
        
        def E_score():
            '''
            Calculates the e_score for the clusterd genes.
            preconditions:  /
            postconditions: returns the e_score as an integer.
            '''
            
            mean = 0

            for k in range(1,self.k):
                score = np.sum([np.square(abs(np.subtract(self.centroids[k], self.lib_norm[dis]))) for dis in self.labeled_data[k]])
                mean += score

            e_score = mean / self.k

            return e_score

        def normalize():
            '''
            Normalises the coordinates of all genes.
            preconditions:  /
            postconditions: sets the normalised coordinates as content for the variable self.lib_norm
            '''
            
            normalized = dict()
            for i in data:
                norm = np.sqrt(np.sum(np.square(np.array(data[i]).astype(float)))) ## SET DATA TYPE TO FLOAT --> THEN CALCULATE THE LENGTH OF THE VECTOR BY APPLYING THE FORMULA IN THE CASUS
                normalized[i] = np.divide(np.array(data[i]).astype(float), norm) ## DIVIDE EACH VALUE IN THE COORDINATES BY THE LENGTH OF THE VECTOR
            
            self.lib_norm = normalized
        
        
        def rand_label():
            '''
            Assigns random genes to clusters.
            preconditions:  /
            postconditions: returs a library:
                                keys   --> cluster
                                values --> genes that are randomly assigned to that cluster
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
            Calculates the centers of the clusters.
            preconditions:  /
            postconditions: returns a library:
                                keys   --> cluster
                                values --> list of coordinats for the center of that cluster
            '''
            
            centroids = dict()
            
            # CALCULATE THE MEAN VECTOR PER CLUSTER AND SAVE IT AS CENTROID OF THE CLUSTER
            for k in range(1,self.k):
                centroids[k] = np.mean(np.array([self.lib_norm[index] for index in self.labeled_data[k]]), axis=0) 

            return centroids
        
        def assign_cluster(unlabeld_data):
            '''
            Assigns genes to the cluster with the nearest center.
            preconditions:  a library:
                                keys  --> gen
                                value --> the normalised coordinates of the gen
            postconditions: a library:
                                keys  --> cluster
                                value --> genes assigned to that cluster                     
            '''
            
            lib_clust = {k:[] for k in range(1,self.k)}
            
            for coords in unlabeld_data:
                # GET THE CLUSTER NUMBER WHERE THE DISTANCE IS THE SMALLEST FROM VECTOR TO CENTROID
                # np.armin --> GIVES THE INDEX OF THE LOWEST VALUE IN A LIST. THIS RANGES FROM 0 to k-1, SO 1 NEEDS TO BE ADDED AT THE END
                clust = np.argmin([np.sqrt(np.sum(np.square(np.subtract(self.centroids[k],unlabeld_data[coords])))) for k in range(1,self.k)])+1
                lib_clust[clust].append(coords)

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
        Clusters the data multiple times with different initial clusterings, and chooses the best result.
        preconditions:  min_seed,max_seed: range of which seeds to try
        postconditions: a library:
                            keys  --> cluster
                            value --> genes assigned to that cluster                     
        '''
        
        # START WITH THE BIGGEST POSSIBLE #NR SO THE FIRST TESTED SEED IS ALWAYS LOWER THAN THIS #NR SO THE IF STATEMENT WORKS CORRECTLY
        old_e = float('inf') 
        
        # TEST THE SEEDS FOR WHICH THE RANGE IS SET:
        for seed in range(min_seed,max_seeds+1):
            self.seed = seed
            self.fit(self.input_data)
            new_e = self.E_score
            
            # COMPARE THE NEWLY GENERATED E-SCORE WITH THE OLD-BEST E-SCORE --> IF NEWLY GENERATED E-SCORE IS BETTER THEN THE LAST ONE, SAVE IT
            if new_e < old_e: 
                best_labeld_data = self.labeled_data
                best_seed = seed
                old_e = new_e
                
        self.seed = best_seed
        
        return best_labeld_data, best_seed, old_e
        
    
    def predict(self, input_data):
        '''
        Assigns new input data to its nearest cluster, according to previously determined centroids that are not affected by this new data.
        preconditions: input_data is either a dictionary of datapoints or a single datapoint (a list of coordinates)
        postconditions: a library:
                             keys  --> cluster
                             value --> genes assigned to that cluster
               
        '''
        lib_clust = {k:[] for k in range(1,self.k)}
        
        # CHECK IF INPUT TYPE CONTAINS MULTIPLE DATA POINTS OR JUST A SINGLE DATA POINT:
        if type(input_data) == dict:
            for coords in input_data:
                # NORMALIZE INPUT DATA
                norm = np.sqrt(np.sum(np.square(np.array(input_data[coords]).astype(float)))) 
                normalized = np.divide(np.array(input_data[coords]).astype(float), norm)

                # CALCULATE CLOSEST CENTROID AND SAVE IT TO CORRESPONDING CLUSTER
                clust = np.argmin([np.sqrt(np.sum(np.square(np.subtract(self.centroids[k],normalized)))) for k in range(1,self.k)])+1
                lib_clust[clust].append(coords)
                
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
    def __init__(self, r=1, j=10, thresh=4):
        '''
        preconditions:  The range of neighbours (r), the total amount of sub-spaces per axis (j), the threshhold on how many 
                        points there need to be per cell to count it as clustercell (thresh).
        postconditions: sets the range to self.r
                        sets the amount of sub-spaces to self.j
                        sets the threshhold to self.thresh
        '''
        self.r = r
        self.j = j #!!! Moet andere variabelnaam krijgen
        self.thresh = thresh    
    
    def fit(self, data):
        '''
        Assigns all genes to grid cells.
        preconditions:  a library with that contains the data. the library must have the gene code as key and a list
                        of the datapoints as value.
        postconditions: returns a library: 
                            keys   --> clusters
                            values --> unique cell ID
        '''
        self.indexes = list(data.keys())

        def get_coordinates(data):
            '''
            Returns a library with as key every axis and as value a tuple of all value of the axis.
            preconditions:  a library with that contains the data. the library must have tehe gene code as key and a list
                            of the datapoints as value.
            postconditions: returns a library: 
                                keys   --> axis
                                values --> list of all values in that specific axis
            '''
            coordinates = dict()

            for cor in range(len(data[self.indexes[0]])): # GET THE GENERAL LENGTH OF THE COORDINATES --> SO FOR EVERY AXIS THERE WILL BE A LIST
                x = list()
                for i in data: # GET EVERY PROTEIN ID IN THE DATA DICTIONARY
                    x.append(float(data[i][cor]))  # APPEND THE AXIS LIST WITH THE CORROSPONDING AXIS FROM EVERY PROTEIN

                coordinates[cor] = tuple(x) # MAKE TUPLE OF THE AXIS VALUES
            return coordinates

        def get_intervals(data,j):
            '''
            Returns a library with the subspaces of the data sorted on the axis.
            preconditions:  a library that contains the axis-data.
            postconditions: return per axis the intervals of the subspaces
            '''
            dic = dict() #!!! Misschien andere naam?

            # GET THE RANGE (MIN-MAX) OF EVERY AXIS --> SO SUB-cellS CAN BE MADE IN THIS RANGE
            for i in data:
                minim = np.min(data[i])
                maxim = np.max(data[i])

                dic[i] = np.mgrid[minim:maxim:complex(0,j)] # J IS THE #NR OF INTERVALLS WE WANT EVERY AXIS

            return dic

        def get_cells():
            '''
            Returns a library with gene ID as key and as value the unique cell ID.
            preconditions:  - a library with that contains the data. the library must have tehe gene code as key and a list
                              of the datapoints as value. 
                            - a library with the intervals
                            
            postconditions: returns a library: 
                                keys   --> gene ID
                                values --> unique cell ID
            '''
            coordinates = get_coordinates(data) # GET AXIS DATA
            intervals = get_intervals(coordinates,self.j) # GET INTERVAL DATA FROM AXIS DATA --> TO MAKE CONDITIONS FOR cell GENERATION
            
            labels = list()
            
            # cell-ID GENERATION VERTICALLY:
            # --> FOR EVERY AXIS A LIST IS GENERATED WITH THE SUBSPACE INDEX (THIS IS A #NR WHICH RANGES FROM 0 TO THE MAX NUMBER OF SUBSPACES THAT CAN BE MADE)
            for i in coordinates: # LOOP OVER EVERY AXIS
                x_inter_label = list() # MAKE NEW LIST FOR EVERY AXIS
                for j in range(1,len(intervals[i])): # LOOP OVER THE CONDITIONS FOR EVERY AXIS --> INTERVALS
                    # MAKE A LIST WHICH CONTAINS OF TRUE OR FALSE STATEMENS --> DUE TO THE CONDITIONS THAT ARE IMPOSED
                    # --> A LIST OF THE LENGTH OF EVERY AXIS IS RETURNED
                    vals = list((coordinates[i] >= intervals[i][j-1]) & (coordinates[i] < intervals[i][j])) # MAKE BOOLEAN LIST 
                    vals = np.where(vals==False, 0, vals) # IF CONDITIONS ARE NOT MET SET VALUE TO 0
                    vals = np.where(vals==True, j, vals) # IF CONDITIONS ARE MET SET VALUE TO THE CORRESPONDING SUBSPACE THE VALUE BELONGS TO      

                    x_inter_label.append(vals) # APPEND THE LIST WITH THE GATHERED SUBSPACE INDEXES

                labels.append(x_inter_label) # APPEND THE LABELS LIST WITH THE LISTS CONTAINTING THE INDEXES OF THE SUBSPACES PER AXIS


            # REWRITE THE LABELS LIST TO A LIST THAT CONTAIS A LIST WITH THE CORRESPONDING SUBSPACE INDEXES PER PROTEIN
            # THE LABELS LIST CONTAINS PER AXIS A NUMBER OF LISTS FOR EVERY SUB-SPACE. TO COMBINE THESE LISTS, THEY CAN SIMPLY BE ADDED TO EACHOTHER
                # IF THE DATA DID NOT SATISFY THE CONDITIONS THE #NR 0 WAS AWARDED 
                # IF THE DATA DID SATISFY THE CONDITIONS THE #NR OF THE SUBSPACE WAS AWARDED
            # --> THE NR OF THE SUBSPACE THAT IS CORRECTLY AWARDED TO THE VALUE IS NOT INFLUENCED WHEN WE ADD 0 TO IT (EVEN IF THE SUBSPACE AWARDED WAS 0, THEN IT REMAINS 0)
            # --> SO IF WE ADD UP THE LIST OF ALL SUBSPACE OPTIONS WE GET A SINGLE LIST OF SUBSPACES PER PROTEIN ID PER AXIS
            labelzz = list()
            for i in range(len(labels)):
                for k in range(1,len(labels[i])):
                    labels[i][0] = np.add(labels[i][0],labels[i][k])
                labelzz.append(labels[i][0])

            del labels

            self.lib_cells = dict()     
            
            # SAVE THE LABELZZ LIST IN A DICTIONARY WITH AS KEY THE PROTEIN ID AND AS VALUE THE cell-ID WICH IS A LIST CONTAINING THE INDEXES OF THE SUBSPACES
            for i in range(len(labelzz[0])):
                cell_id = list()
                for j in range(len(labelzz)):
                    cell_id.append(labelzz[j][i])

                self.lib_cells[self.indexes[i]] = cell_id
            
            
            return self.lib_cells
        return get_cells()
        
    def self_predict(self, cell_ids):
        '''
        Assign all genes to their right clusters.
        preconditions:  a library with that contains the data. the library must have the gene code as key and a list
                        of the datapoints as value.
        postconditions: returns a library: 
                            keys   --> clusters
                            values --> genes that are assigned to that cluster
        '''       
        new_indexes = self.indexes.copy()

        def get_neighbour_ID(ID):
            '''
            Calculate which indexes are in neighbour cells of the given cell.
            preconditions: 'ID' is the id of a gene
            postconditions: returns a list of the id's of the genes that are neighbours of the input gene
            '''
            neighbours = list()
            for index in new_indexes:
                distance = np.subtract(self.lib_cells[ID],self.lib_cells[index],out=np.array([0 for i in range(len(self.lib_cells[ID]))])) # CALCULATE THE DISTANCE OF A POTENTIAL NEIGHBOUR TO THE CURRENT CELL
                if all([-self.r <= i <= self.r for i in distance]): # CHECK IF ALL THE VALUES (= ALL DIMENSIONS) OF THE SUBSTRACTION ARE WITHIN THE RANGE
                    neighbours.append(index)
            return neighbours

        def group_neighbours(ID,main_cluster):
            '''
            Places a gene and each of its neighbours in the same cluster.
            precondition: 'ID' is the id of a gene
                          'main_cluster' is a list
            postcondition: appends every gene that should be placed in the main cluster to the list 'main_cluster'
            '''
            neighbours = get_neighbour_ID(ID)
            for neighbour in neighbours: # FIRSTLY, EVERY NEIGHBOUR SHOULD BE ADDED TO THE CLUSTER AND REMOVED FROM new_indexes TO AVOID UNNECESSARY REPETITION (THE SAME NEIGHBOURS BEING DETECTED MULTIPLE TIMES)
                main_cluster.append(neighbour)
                new_indexes.remove(neighbour)
            if len(neighbours) > 1: # IF A GENE HAS EXACTLY ONE NEIGHBOUR, IT IS THE GENE ITSELF. IT IS UNNECESSARY TO CALL group_neighbours() IN THAT CASE.
                 for neighbour in neighbours:
                     group_neighbours(neighbour,main_cluster) # THIS FUNCTION IS RECURSIVE (= CALLS ITSELF); THE PROCESS OF FINDING NEIGHBOURS AND PLACING THEM IN THE MAIN CLUSTER IS REPEATED FOR NEIGHBOURS OF NEIGHBOURS ...ETC.
            return
        
        self.clustered_data = dict()
        cluster_nr = 0
        anticluster_nr = -1
        for index in list(new_indexes): # BECAUSE THIS LOOP REMOVES VALUES FROM new_indexes, A COPY OF THE LIST MUST BE MADE TO ITERATE OVER ALL VALUES. OTHERWISE, SOME INDEXES DISAPPEAR.
            if index in new_indexes:
                cluster = []
                group_neighbours(index,cluster)
                if len(cluster) > self.thresh:                  # CHECK IF GROUPS CONTAIN ENOUGH DATAPOINTS TO BE APPOINTED TO A CLUSTER
                    self.clustered_data[cluster_nr] = cluster   #
                    cluster_nr += 1                             #
                else:                                           # IF NOT --> GIVE IT AN ANTI-CLUSTER --> SO WE CAN DISTINGUISH THE DATA THAT IS NOT ASSIGNED TO A CLUSTER/GROUP 
                    self.clustered_data[anticluster_nr] = cluster     
                    anticluster_nr -= 1
        return self.clustered_data


"""

    GET RESULTS BY CALLING CLASSES
    
"""
        
lib_data = file_to_lib('Data\Voorbeeld_clusterdata.txt')
lib_results = file_to_lib('Data\Voorbeeld_clusterresult.txt')

kmca = KMCA()
kmca.fit(lib_data)
kmca_result, kmca_best_seed, kmca_e_score = kmca.optimize(10,30)

gca = GCA(j=40)
gca_fit_data = gca.fit(lib_data)
gca_results = gca.self_predict(gca_fit_data)


"""

    PRINT TEXT FILES FORMATTED/SORTED LIKE THE INPUT FILE
    
"""

def return_txt_file(data, name, format_data=lib_data):
   '''
   Converts the input data to a text file.
   preconditions:  data is a library with all clusters as keys, and the genes assigned to each key as values;
                   name is a string of the desired name of the text file;
                   format_data library, indicating how to format the data to a text file
   postconditions: returns a text file
   '''
   lib_unfolded = dict()
    
   for cluster in data:
       for ID in data[cluster]:
           lib_unfolded[ID] = cluster
    
   #print(lib_unfolded)
   file = open(f'{name}results.txt', 'w')
   indices_lost = list()
    
   for INDEX in format_data:
       try:
           line = str(INDEX) + " " + str(lib_unfolded[INDEX]) + "\n"
           file.write(line)
       except:
           indices_lost.append(INDEX)
    
   print(f'Function: {name} --> these indices are lost: {indices_lost}')
        
   file.close()
    
#return_txt_file(kmca_result, 'kmca_')
return_txt_file(gca_results, 'gca_')