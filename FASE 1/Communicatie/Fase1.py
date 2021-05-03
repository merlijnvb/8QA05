#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 29 10:17:49 2021

@author: irisalmekinders
"""

#libraries
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import math

#global variables
Columns = ['P1Sig', 'P1STB', 'P1Cov', 'P2Sig', 'P2STB', 'P2Cov', 'P2SigNorm', 'Log_P1Sig', 'Log_P2Sig', 'Log_P2SigNorm']
Day_numbers = [1,2,4,7,14,21,45,90]

def readfile(filename):
    """Preconditions:  Filename is de naam van een .txt bestand dat te vinden
                        is in dezelfde map als deze functie.
    Postconditions:   Retourneert een dictionary met daarin per ID een lijst 
                        van alle bijbehorende waarden."""
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
    """Preconditions:  Alle aan te roepen .txt bestanden staan in dezelfde
                        folder als deze functie.
    Postconditions:   Retourneert een lijst van dictionaries van alle dagen"""
    return [readfile("dag1.txt"),readfile("dag2.txt"),readfile("dag4.txt"),readfile("dag7.txt"),readfile("dag14.txt"),readfile("dag21.txt"),readfile("dag45.txt"),readfile("dag90.txt")]


def summation(Dag_dict,column_head):
    """Preconditions:  Dag_dict is een dictionary met per ID de bijbehorende 
                        waarden en column head is de index van de kolom die 
                        moet worden opgesomd.
    Postconditions: Retourneert de som van alle waarden uit de kolom."""
    total = 0
    
    for key in Dag_dict:
        total += int(Dag_dict[key][Columns.index(column_head)])
    return total

def logaritm(Dag_dict):
    """Preconditions:  Dag_dict is een dictionary met per ID de bijbehorende waarden.
    Postconditions:  Zet P1Sig, P2Sig en P2SigNorm in een logaritme en voegt deze toe aan de dictionary."""
    for key in Dag_dict:
        Dag_dict[key].append(math.log10(Dag_dict[key][0]))
        Dag_dict[key].append(math.log10(Dag_dict[key][3]))
        Dag_dict[key].append(math.log10(Dag_dict[key][6]))
        
def Log_add(Dagen):
    for Dag in Dagen:
        logaritm(Dag) # logaritm functie aanroepen, zodat de logaritmische waardes worden toegevoegd
   
def normalize(Dag_dict):
    """Preconditions:  Dag_dict is een dictionary met per ID de bijbehorende waarden.
    Postconditions:  Normaliseert P2Sig (P2Sig * S1/S2) en voegt deze toe aan de dictionary."""
    S1 = summation(Dag_dict, "P1Sig") # optellen van alle P1Sig waardes
    S2 = summation(Dag_dict, "P2Sig") # optellen van alle P2Sig waardes
    for key in Dag_dict:
        Dag_dict[key].append(int(Dag_dict[key][3])*S1/S2)

def normSig2_add(Dagen):
    for Dag in Dagen:
        normalize(Dag) # normalize functie aanroepen, zodat de normalisatie waardes worden toegevoegd
        
def plot_dag(ax_dag,df_dag,x,y,log=False):
    """Preconditions: dict van een dag, plotnummer, log=T/F,x/y kolommen weten"""
    df_dag.plot(kind='scatter', x=x, y=y, c='r', ax=ax_dag)
    

def plot_phase1(Dagen,log=False,Norm=False):
    #make df
    add_on = ''
    norm_fac = 0
    if log:
        add_on += ', logaritmisch'
        if Norm:
            add_on += ' en genormaliseerd'
            norm_fac = 1  
    elif Norm: 
        add_on += ', genormaliseerd'
        norm_fac = 3
    
    fig, ax = plt.subplots(2,4,figsize=(20,10),sharex=True,sharey=True)
    fig.suptitle("Visualisatie ruwe data"+add_on,size=24,weight='bold')
    for i in range(len(Dagen)):
        ax_dag = ax[i%2,i//2]
        ax_dag.set_title("Dag "+str(Day_numbers[i]))
        df_dag = pd.DataFrame.from_dict(Dagen[i]).transpose()
        df_dag.columns = Columns
        if log:
            plot_dag(ax_dag,df_dag,Columns[7],Columns[8+norm_fac], True)
            ax_dag.plot([0,5],[0,5], c="k")
        else: 
            ax_dag.plot([0,50000],[0,50000], c="k")
            plot_dag(ax_dag,df_dag,Columns[0],Columns[3 + norm_fac],True)


def Main():
    Dagen = make_Days()
    normSig2_add(Dagen)
    Log_add(Dagen)
    plot_phase1(Dagen)
    plot_phase1(Dagen,log=True)
    plot_phase1(Dagen, Norm = True)
    plot_phase1(Dagen, log = True, Norm = True)

Main()


