import pandas as pd
import seaborn as sns

def lib(FileName):
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
            
lib_data = lib('Voorbeeld_clusterdata.txt')
lib_results = lib('Voorbeeld_clusterresult.txt')

def plot_clusters(data, results):
    df_data = pd.DataFrame(data)
    df_data = df_data.transpose()
    clusters = list(results.values())
    df_data['cluster'] = clusters
    
    for cluster in range(df_data['cluster'].max()+1):
        df_cluster = df_data[df_data['cluster']==cluster]
        df_cluster = df_cluster.drop('cluster', axis=1)
        df_cluster = df_cluster.astype(float)
        df_cluster = df_cluster.transpose()
        #ax = df_cluster.plot(legend=False)
        ax = df_cluster.mean(axis=1).plot(legend=False)
        #ax.set_title(cluster)

plot_clusters(lib_data, lib_results)

    


