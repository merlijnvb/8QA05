# -*- coding: utf-8 -*-
"""
Created on Thu Apr 29 10:10:55 2021

@author: 20202309
"""

#libraries
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

#global variables
Columns = 'P1Sig	P1STB	P1Cov	P2Sig	P2STB	P2Cov\tP2SigNorm'.split('\t')

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
        
def plot_dag(df_dag):
    'hello'
    

def plot_phase1(Dagen):
    #make df
    for Dag in Dagen[:2]:
        df_dag = pd.DataFrame.from_dict(Dag).transpose()
        df_dag.columns = Columns
        plot_dag(df_dag)
    
    

def Main():
    Dagen = make_Days()
    normSig2_add(Dagen)
    plot_phase1(Dagen)
    


Main()

