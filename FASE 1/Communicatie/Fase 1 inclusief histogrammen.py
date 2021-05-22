#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Thu Apr 29 10:17:49 2021
@author: Iris Almekinders en Niels van Noort
"""

#libraries
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import math

#global variables
Columns = ['P2SigNorm', 'Log_P1Sig', 'Log_P2Sig', 'Log_P2SigNorm', 'Class']
Day_numbers = []

def find_files(file_key_word):
    """returns list of textfilenames based on Day_numbers"""
    from os import listdir
    from pathlib import Path
    
    arr = listdir(Path(__file__).parent.absolute())
    filenames = [filename for filename in arr if f'{file_key_word}' in filename]
    
    day_nrs = [name.replace('.txt','') for name in filenames]
    day_nrs = [int(name.replace(f'{file_key_word}','')) for name in day_nrs]
    day_nrs.sort()
    
    filenames = []
    for number in day_nrs: 
        Day_numbers.append(number)
        filenames.append("dag"+str(number)+".txt")
        
    return [readfile(filename) for filename in filenames]        

def readfile(filename):
    """Preconditions:  Filename is de naam van een .txt bestand dat te vinden
                        is in dezelfde map als deze functie.
    Postconditions:    Retourneert een dictionary met daarin per ID een lijst 
                        van alle bijbehorende waarden."""
    #read the text file
    infile = open(filename)
    lines = infile.readlines()
    infile.close()
    ID_dict = {}
    if len(Columns) == 5: 
        for col_head in reversed(lines.pop(0).split()[1:]): Columns.insert(0,col_head) # adds the headings to global var Columns, only if that hasn't been done yet
    else: lines.pop(0)
    
    # make a dictionary-entry per line
    for line in lines:
        values = line.rstrip().split('\t')
        
        # make all values integers
        str_vals = values[1:]
        for i in range(len(str_vals)):
            str_vals[i] = int(str_vals[i])
            
        # add line to list
        ID_dict[int(values[0])] = str_vals
    return ID_dict

def summation(Dag_dict,column_head):
    """Preconditions:  Dag_dict is een dictionary met per ID de bijbehorende 
                        waarden en column head is de index van de kolom die 
                        moet worden opgesomd.
    Postconditions:    Retourneert de som van alle waarden uit de kolom."""
    total = 0 # counter
    
    for key in Dag_dict:
        total += int(Dag_dict[key][Columns.index(column_head)])
    return total

def Log_add(Dagen):
    """Preconditions:  dagen is een lijst van dictionaries met per ID de 
                        bijbehorende waarden.
    Postconditions:    Zet P1Sig, P2Sig en P2SigNorm in een logaritme en voegt 
                        deze toe aan elke dictionary in Dagen."""
    for Dag_dict in Dagen:
        for key in Dag_dict:
            Dag_dict[key].append(math.log10(Dag_dict[key][0]))
            Dag_dict[key].append(math.log10(Dag_dict[key][3]))
            Dag_dict[key].append(math.log10(Dag_dict[key][6]))

def normSig2_add(Dagen):
    """Preconditions:  Dagen is een lijst van dictionaries met per ID de 
                        bijbehorende waarden.
    Postconditions:    Normaliseert P2Sig (P2Sig * S1/S2) en voegt deze toe aan
                        elke dictionary."""
    for Dag in Dagen:
        S1 = summation(Dag, "P1Sig") # optellen van alle P1Sig waardes
        S2 = summation(Dag, "P2Sig") # optellen van alle P2Sig waardes
        for key in Dag:
            Dag[key].append(int(Dag[key][3])*S1/S2)

def plot_dagen(Dagen,log=False,Norm=False):
    """Preconditions:  Dagen is een lijst van libraries, log en
                        norm zijn booleans die aangeven of de 
                        functie logaritmisch, dan wel genormaliseerd
                        moet worden geplot.
    Postconditions:    Zorgt ervoor dat de data van elke dag in één
                        figuur zou """    
    # make extra variables
    add_on = '' # for the title
    norm_fac = 0 # for normalisation column index
    Colors = {'A': 'red', 'B': 'blue', 'C': 'green', 'D': 'yellow'}
    
    # use the booleans for the title
    if log:
        add_on += ', logaritmisch'
        if Norm:
            add_on += ' en genormaliseerd'
            norm_fac = 1  
    elif Norm: 
        add_on += ', genormaliseerd'
        norm_fac = 3
    
    # start plotting
    fig, ax = plt.subplots(2,4,figsize=(20,10),sharex=True,sharey=True)
    fig.suptitle("Visualisatie ruwe data"+add_on,size=24,weight='bold')
    for i in range(len(Dagen)):
        ax_dag = ax[i%2,i//2]
        ax_dag.set_title("Dag "+str(Day_numbers[i]))
        df_dag = pd.DataFrame.from_dict(Dagen[i]).transpose()
        df_dag.columns = Columns
        # assign color to right class for plotting
        df_dag = df_dag.replace({'Class': Colors})
        # plot every class per day; this ensures the legenda can be made
        for Class in Colors:
            # making a DataFrame for every class using a mask
            df_dag_new = df_dag[df_dag['Class'] == Colors[Class]]
            # check whether we should use the logaritmic columns
            if log:
                # plotting one class
                df_dag_new.plot(kind='scatter', x=Columns[7], y=Columns[8+norm_fac], c=Colors[Class], ax=ax_dag, colorbar = False)
                # only plot P1 = P2 once, not four times
                if Class == 'A': ax_dag.plot([1,5],[1,5], c="k")
            else: 
                # only plot P1 = P2 once, not four times
                if Class == 'A': ax_dag.plot([0,50000],[0,50000], c="k")
                # plotting one class
                df_dag_new.plot(kind='scatter', x=Columns[0], y=Columns[3 + norm_fac], c=Colors[Class], ax=ax_dag, colorbar = False)
            ax_dag.legend(['P1=P2',"A","B","C","D"])
     
def plot_phase1(Dagen):
    """Preconditions:  Dagen is een lijst van libraries
    Postconditions:    Roept alle mogelijkheden voor  
                        plot_dagen aan."""   
    for log in [False,True]:
        for norm in [False,True]:
            plot_dagen(Dagen,log,norm)
            
def make_class(Dag,stb_threshold=25):
    """Preconditions:  Dag is een dictionary met per ID de bijbehorende 
                        waarden. stb_threshold is de drempelwaarde voor de 
                        filtering van de classes.
    Postconditions:    Voegt aan elke waarde van Dag een klasse (A, B, C of D)
                        toe op basis van de waarden van P1STB en P2STB. Voor
                        deze klassen geldt respectievelijk in vergelijking
                        met de drempelwaarde: Geen waarden lager, slechts P2 
                        lager, slechts P1 lager, beide lager.
                        Retourneert een lijst met de frequenties van alle
                        klassen"""
    counts = [0,0,0,0]
    
    for ID in Dag:
        if (Dag[ID][1] >= stb_threshold):
            if (Dag[ID][4] >= stb_threshold):
                Dag[ID].append("A")
                counts[0]+=1
            else:
                Dag[ID].append("B")
                counts[1]+=1
        else:
            if (Dag[ID][4] >= stb_threshold):
                Dag[ID].append("C")
                counts[2]+=1
            else:
                Dag[ID].append("D")
                counts[3]+=1
    return counts   
            
def Classes(Dagen,frequences=False):
    """Preconditions:  Dagen is een lijst van libraries, frequences is een 
                        boolean die het printen aan kan zetten, per default uit.
    Postconditions:    Roept de funtie make_class aan per dag en print daarmee
                        de frequenties per dag en klasse als frequences==True"""
    classes = "ABCD"
    
    for i in range(len(Dagen)):
        counts = make_class(Dagen[i])
        if frequences: 
            print("\nVoor dag", Day_numbers[i], "geldt: ")
            for j in range(len(counts)):
                print("\tFrequentie van", classes[j],"is",counts[j]) 


def add_expression(Dagen):
    """Voegt de relatieve expressiewaarde toe aan de bestaande dictionaries."""
    Columns.append("RelativeExpression")
    for Dag in Dagen:
        for ID in Dag:
            if (Dag[ID][0]/Dag[ID][6]) >= 1:
                r = (Dag[ID][0]/Dag[ID][6]) - 1
            else:
                r = (-Dag[ID][6]/Dag[ID][0]) + 1
            Dag[ID].append(r)

def Daysdict(Dagen,r_filter=0.5):
    """New_Dagen is een dictionary met per ID de r voor elke dag"""
    New_Dagen = {}
    Filterinfo = {}
    
    for ID in Dagen[0]:
        Rvals = list()
        filterinfo = []
        Filterinfo[ID] = [0,0,0,0]
        
        for i in range(len(Day_numbers)):
            Rvals.append(Dagen[i][ID][-1]) 
            if -r_filter < Dagen[i][ID][-1] < r_filter: Filterinfo[ID][3] +=1 # counter voor hoevaak r < filterwaarde
            filterinfo.append(Dagen[i][ID][10])
            if not ((Dagen[i][ID][2] >= 40) and (Dagen[i][ID][2] <= 160)): Filterinfo[ID][2] += 1
    
        for class_name in filterinfo:
            if (class_name == "B") or (class_name == "C"): Filterinfo[ID][0] += 1
            if class_name == "D": Filterinfo[ID][1] += 1
            
        New_Dagen[ID] = Rvals
    return New_Dagen, Filterinfo
    
def bar_plots(rDict,r_filter=0.5,movement = 256):
    """ grootste voortgang van ZSO 5
    maakt 256 hele mooie barplots --> nutteloos voor fase 1 :) """
    # ncols = math.ceil(len(rDict)**0.5)
    ncols = 10
    nrows = ncols
    fig, ax = plt.subplots(ncols,nrows,figsize=(40,40),squeeze=False,sharex=True,sharey=True)
    fig.suptitle("Fucking veel lijndiagrammen",size=60,weight='bold') # this will be changed... for now we're a bit pissed at the number of plots and the time it takes to draw them
    df_expr = pd.DataFrame.from_dict(rDict)
    for col in range(ncols):
        for row in range(nrows):
            ID = df_expr.columns[col*ncols+row+movement]
            df_expr[ID].plot(ax=ax[row,col])
            ax[row,col].set_title(str(ID))
            ax[row,col].axhline(y=r_filter)
            ax[row,col].axhline(y=-r_filter)

def filtering(rDict, Filterinfo):
    """er komt Dagen in en gaat vervolgens per ID kijken hoevaak classes B en c
    voorkomen, hoevaak de spotgrootte te laag is en of r ooit boven de r_filter
    uitkomt."""
    filtered_rDict = {}
    
    for ID in rDict:
        if (Filterinfo[ID][0] + Filterinfo[ID][2] < 3) and (Filterinfo[ID][1] == 0) and (Filterinfo[ID][3] < 8):
            filtered_rDict[ID] = rDict[ID]
    return filtered_rDict

def file_phase3(rDict):
    outf = open('Ruwe Data Fase 3.txt', 'w')
    
    for ID in rDict:
        outf.write(str(ID))
        outf.write('\t')
        r_values = []
        for i in rDict[ID]:
            r_values.append(str(i))
        outf.write('\t'.join(r_values))
        outf.write('\n')
        
    outf.close()
    
def plot_hist(rDict):
    """""Plot een histogram van de expressiewaarden voor elke dag, zodat er
    te zien is hoe vaak een expressiewaarde voorkomt."""
    fig, ax = plt.subplots(2,4,figsize=(20,10),sharex=True,sharey=True)
    fig.suptitle("Histogram expressiewaarden per dag", size=24,weight='bold')
    
    df = pd.DataFrame.from_dict(rDict).transpose()

    for i in range(len(Day_numbers)):
        ax_dag = ax[i%2,i//2]
        ax_dag.set_title("Dag "+str(Day_numbers[i]))
        ax_dag.set_xlim(-5, 5)
        df[i].plot(kind = 'hist', ax = ax_dag)
        
    
def Main():
    Dagen = find_files('dag')
    normSig2_add(Dagen)
    Log_add(Dagen)
    Classes(Dagen,frequences=True)
    add_expression(Dagen)
    rDict, Filterinfo = Daysdict(Dagen)
    filtered_rDict = filtering(rDict, Filterinfo)
    #file_phase3(rDict)
    plot_hist(rDict)
    # bar_plots(rDict,0.5)
    # plot_phase1(Dagen)
    
Main()




