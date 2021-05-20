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
                normalized[i] = np.divide(np.array(data[i]).astype(float), norm) ## DIVIDE EACH VALUE IN THE CORDINATES BY THE LENGTH OF THE VECTOR
            
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
        Assigns all genes to grid cells.
        preconditions:  a library with that contains the data. the library must have the gene code as key and a list
                        of the datapoints as value.
        postconditions: returns a library: 
                            keys   --> clusters
                            values --> unique cell ID
        '''
        self.indexes = list(data.keys())

        def get_cordinates(data):
            '''
            Returns a library with as key every axis and as value a tuple of all value of the axis.
            preconditions:  a library with that contains the data. the library must have tehe gene code as key and a list
                            of the datapoints as value.
            postconditions: returns a library: 
                                keys   --> axis
                                values --> list of all values in that specific axis
            '''
            cordinates = dict()

            for cor in range(len(data[self.indexes[0]])): # GET THE GENERAL LENGTH OF THE CORDINATES --> SO FOR EVERY AXIS THERE WILL BE A LIST
                x = list()
                for i in data: # GET EVERY PROTEIN ID IN THE DATA DICTIONARY
                    x.append(float(data[i][cor]))  # APPEND THE AXIS LIST WITH THE CORROSPONDING AXIS FROM EVERY PROTEIN

                cordinates[cor] = tuple(x) # MAKE TUPLE OF THE AXIS VALUES
            return cordinates

        def get_intervals(data,j):
            '''
            Returns a library with the subspaces of the data sorted on the axis.
            preconditions:  a library that contains the axis-data.
            postconditions: return per axis the intervals of the subspaces
            '''
            dic = dict()

            # GET THE RANGE (MIN-MAX) OF EVERY AXIS --> SO SUB-CELLS CAN BE MADE IN THIS RANGE
            for i in data:
                minim = np.min(data[i])
                maxim = np.max(data[i])

                dic[i] = np.mgrid[minim:maxim:complex(0,j)] # J IS THE #NR OF INTERVALLS WE WANT EVERY AXIS

            return dic

        def get_cells():
            '''
            Returns a library with gene ID as key and as value the unique cel ID.
            preconditions:  - a library with that contains the data. the library must have tehe gene code as key and a list
                              of the datapoints as value. 
                            - a library with the intervals
                            
            postconditions: returns a library: 
                                keys   --> gene ID
                                values --> unique cel ID
            '''
            cordinates = get_cordinates(data) # GET AXIS DATA
            intervals = get_intervals(cordinates,self.j) # GET INTERVAL DATA FROM AXIS DATA --> TO MAKE CONDITIONS FOR CELL GENERATION
            
            labels = list()
            
            # CELL-ID GENERATION VERTICALLY:
            # --> FOR EVERY AXIS A LIST IS GENERATED WITH THE SUBSPACE INDEX (THIS IS A #NR WHICH RANGES FROM 0 TO THE MAX NUMBER OF SUBSPACES THAT CAN BE MADE)
            for i in cordinates: # LOOP OVER EVERY AXIS
                x_inter_label = list() # MAKE NEW LIST FOR EVERY AXIS
                for j in range(1,len(intervals[i])): # LOOP OVER THE CONDITIONS FOR EVERY AXIS --> INTERVALS
                    # MAKE A LIST WHICH CONTAINS OF TRUE OR FALSE STATEMENS --> DUE TO THE CONDITIONS THAT ARE IMPOSED
                    # --> A LIST OF THE LENGTH OF EVERY AXIS IS RETURNED
                    vals = list((cordinates[i] >= intervals[i][j-1]) & (cordinates[i] < intervals[i][j])) # MAKE BOOLEAN LIST 
                    vals = np.where(vals==False, 0, vals) # IF CONDITIONS ARE NOT MET SET VALUE TO 0
                    vals = np.where(vals==True, j, vals) # IF CONDITIONS ARE MET SET VALUE TO THE CORRESPONDING SUBSPACE THE VALUE BELONGS TO      

                    x_inter_label.append(vals) # APPEND THE LIST WITH THE GATHERD SUBSPACE INDEXES

                labels.append(x_inter_label) # APPEND THE LABELS LIST WITH THE LISTS CONTAINTING THE IDNEXES OF THE SUBSPACES PER AXIS


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
            
            # SAVE THE LABELZZ LIST IN A DICTIONARY WITH AS KEY THE PROTEIN ID AND AS VALUE THE CELL-ID WICH IS A LIST CONTAINING THE INDEXES OF THE SUBSPACES
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
        def get_neigb_id(ID):
            '''
            Returns neighbour cell ID's.
            preconditions:  a cell ID
            postconditions: returns a list with all the IDs of the neigboring cells
            '''
            neigb = list()
            #GET EVERY NEIGBOUR OF THE INPUT-CELL:
            for indx in cell_ids: # LOOP OVER EVERY CELL-INDEX THAT IS GENERATED
                # THE DISTANCE TO THE INPUT-CELL IS A FORM OF MANHATTAN DISTANCE --> SUBTRACT EVERY SUBSPACE-ID FROM THE INPUT-CELL SUBSPACES
                # --> WHEN THE OUTCOME OF EVERY SUBSPACE OF A CELL-ID IS IN THE SET RANGE, IT CAN BE SEEN AS A NEIGBOR
                distance = np.subtract(self.lib_cells[ID],self.lib_cells[indx],out=np.array([0 for i in range(8)])) 
                if all([-self.r < i < self.r for i in distance]): # CHECK IF ALL THE VALUES OF THE SUBSTRACTION (in distance) ARE IN RANGE
                    neigb.append(indx) # APPEND THE NEIGBOUR LIST WITH THE PROTEIN ID OF THE NEIGBOUR

            return neigb

        new_indexes = self.indexes.copy() # MAKE DUMMY LIST OF THE PROTEIN IDs

        def make_groups(data_point):
            '''
            Returns a group/cluster.
            preconditions:  a cell ID
            postconditions: returns a list with all the neigbours of neigbours that originated from the first cell
            '''
            neigbors = get_neigb_id(data_point) # GET THE FIRST NEIGBORS OF A DATAPOINT (PROTEIN ID)
            new_neigbors = list()

            for cel in neigbors: # GATHER FOR EVERY NEIGBOUR IF IT HAS NEW UNIQUE NEIGBOURS
                new_cells = get_neigb_id(cel)

                for new_cel in new_cells:
                    if new_cel not in neigbors: # CHECK IF NEIGBOUR IS UNIQUE
                        new_neigbors.append(new_cel)
                
                # ADD NEWLY FOUND NEIGBOURS TO THE ORIGINAL NEIGBOURS LIST
                neigbors += new_neigbors

            for neigb in neigbors: 
                # REMOVE ALL THE NEIGBOURS WHICH ARE ALREADY ASSIGNED FROM THE DUMMY PROTEIN INDEX LIST
                new_indexes.remove(neigb)

            return neigbors

        options = list()

        # RUN THE MAKE GROUPS ON THE DUMMY LIST --> SO THAT WE DO NOT HAVE DUPLICATE GROUPS
        # --> MAKES THE PROGRAM FASTER --> NO UNNECESSARY CALCULATIONS
        for index in new_indexes:
            options.append(make_groups(index))


        self.labeled_data = dict()
        cluster = 0
        anti_cluster = -1

        # CHECK IF GROUPS CONTAIN ENOUGH DATAPOINTS TO BE APPOINTED TO A CLUSTER
        # IF NOT --> GIVE IT AN ANTI-CLUSTER --> SO WE CAN DISTINGUISH THE DATA THAT IS NOT ASSIGNED TO A CLUSTER/GROUP 
        for i in range(len(options)):
            if len(options[i]) > self.tresh:
                self.labeled_data[cluster] = options[i]
                cluster += 1
            else:
                self.labeled_data[anti_cluster] = options[i]
                anti_cluster -= 1

        return self.labeled_data



"""

    GET RESULTS BY CALLING CLASSES
    
"""
        
lib_data = file_to_lib('Data\Voorbeeld_clusterdata.txt')
lib_results = file_to_lib('Data\Voorbeeld_clusterresult.txt')

kmca = KMCA()
kmca.fit(lib_data)
kmca_result, kmca_best_seed, kmca_e_score = kmca.optimize(10,30)

gca = GCA()
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
    
return_txt_file(kmca_result, 'kmca_')
return_txt_file(gca_results, 'gca_')
