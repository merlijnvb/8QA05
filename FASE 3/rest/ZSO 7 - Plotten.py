import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import statsmodels.api as sm
from scipy import stats
from scipy.stats import *
import sqlite3
import pandas as pd
import numpy as np


sns.set()
plt.rcParams['figure.figsize'] = 10, 5  # default hor./vert. size of plots, in inches
plt.rcParams['lines.markeredgewidth'] = 1  # to fix issue with seaborn box plots; needed after import seaborn

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

df_data = pd.DataFrame(lib_data)
df_data = df_data.transpose()

clusters = list(lib_results.values())

df_data['cluster'] = clusters
df_data = pd.DataFrame(df_data)
#df_clusters = pd.DataFrame(df_data.groupby('cluster'))
# for cluster in list(df_clusters.index.values):
    

