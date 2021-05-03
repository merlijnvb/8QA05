# -*- coding: utf-8 -*-
"""
Created on Thu Apr 29 10:10:55 2021

@author: 20202309
"""

#libraries
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import math

#global variables
Columns = 'P1Sig	P1STB	P1Cov	P2Sig	P2STB	P2Cov\tP2SigNorm'.split('\t')
Day_numbers = [1,2,4,7,14,21,45,90]

def readfile(filename):
    infile = open(filename)
    lines = infile.readlines()
    infile.close()
    ID_dict = {}
    lines.pop(0)
    
    for line in lines:
        values = line.rstrip().split('\t')
        
        #make all values integers
        str_vals = values[1:]
        for i in range(len(str_vals)):
            str_vals[i] = int(str_vals[i])
            
        #add to list
        ID_dict[int(values[0])] = str_vals
    return ID_dict
    
def make_Days():
    return [readfile("dag1.txt"),readfile("dag2.txt"),readfile("dag4.txt"),readfile("dag7.txt"),readfile("dag14.txt"),readfile("dag21.txt"),readfile("dag45.txt"),readfile("dag90.txt")]


def summation(Dag_dict,column_head):
    total = 0
    
    for key in Dag_dict:
        total += int(Dag_dict[key][Columns.index(column_head)])
    return total
   
def normalize(Dag_dict):
    S1 = summation(Dag_dict, "P1Sig")
    S2 = summation(Dag_dict, "P2Sig")
    for key in Dag_dict:
        Dag_dict[key].append(int(Dag_dict[key][3])*S1/S2)

def normSig2_add(Dagen):
    for Dag in Dagen:
        normalize(Dag)
        
def plot_dag(ax_dag,df_dag,x,y,log=False):
    """Preconditions: dict van een dag, plotnummer, log=T/F,x/y kolommen weten"""
    df_dag.plot(kind='scatter', x=x, y=y, c='r', ax=ax_dag)
    

def plot_phase1(Dagen,log=False,Norm=False):
    #make df
    add_on = ''
    if log:
        add_on += ', logaritmisch'
        if Norm: add_on += ' en genormaliseerd'
    elif Norm: 
        add_on+=', genormaliseerd'
    
    fig, ax = plt.subplots(2,4,figsize=(20,10),sharex=True,sharey=True)
    fig.suptitle("Visualisatie ruwe data I"+add_on,size=24,weight='bold')
    for i in range(len(Dagen)):
        ax_dag = ax[i%2,i//2]
        ax_dag.set_title("Dag "+str(Day_numbers[i]))
        df_dag = pd.DataFrame.from_dict(Dagen[i]).transpose()
        df_dag.columns = Columns
        if log:
            df_dag['Log_P1Sig'] = np.log10(df_dag['P1Sig'])
            df_dag['Log_P2Sig'] = np.log10(df_dag['P2Sig'])
            plot_dag(ax_dag,df_dag,'Log_P1Sig','Log_P2Sig',True)
            ax_dag.plot([0,5],[0,5], c="k")
        else: 
            ax_dag.plot([0,50000],[0,50000], c="k")
            plot_dag(ax_dag,df_dag,Columns[0],Columns[3],True)


def Main():
    Dagen = make_Days()
    normSig2_add(Dagen)
    plot_phase1(Dagen)
    plot_phase1(Dagen,log=True)


Main()

