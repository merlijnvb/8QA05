"""
Course code: 8QA05
group number: 9
members: Iris Almekinders and Niels van Noort
"""


# libraries
import matplotlib.pyplot as plt
import pandas as pd
import math


# global variables
Columns = ['P2SigNorm', 'Log_P1Sig', 'Log_P2Sig', 'Log_P2SigNorm', 'Class']
Day_numbers = []

# definitions
def find_files(file_key_word):
    """Preconditions:  file_key_word is een string die voor kan komen als een
                        deel van een bestandsnaam. Er is minimaal één bestand
                        aanwezig in dezelfde map als dit .py bestand wiens 
                        naam begin met file_key_word, gevolgd door een nummer,
                        waarna direct eindigend met .txt
    Postconditions:    Returns list of textfilenames based on Day_numbers"""
    # importing useful ibraries
    from os import listdir
    from pathlib import Path
    
    # creating a list of filenames
    arr = listdir(Path(__file__).parent.absolute())
    filenames = [filename for filename in arr if f'{file_key_word}' in filename]
    
    # getting the number from the list to sort them
    day_nrs = [name.replace('.txt','') for name in filenames]
    day_nrs = [int(name.replace(file_key_word,'')) for name in day_nrs]
    day_nrs.sort()
    
    # making the final, sorted list
    filenames = []
    for number in day_nrs:
        Day_numbers.append(number)
        filenames.append(file_key_word+str(number)+".txt")
        
    return [readfile(filename) for filename in filenames]        


def readfile(filename):
    """Preconditions:  Filename is de naam van een .txt bestand dat te vinden
                        is in dezelfde map als deze functie.
    Postconditions:    Retourneert een dictionary met daarin per ID een lijst 
                        van alle bijbehorende waarden."""
    # read the text file
    infile = open(filename)
    lines = infile.readlines()
    infile.close()
    ID_dict = {}
    if len(Columns) == 5: 
        for col_head in reversed(lines.pop(0).split()[1:]): Columns.insert(0,col_head) # adds the headings to global var Columns, only if that hasn't been done yet
    else: lines.pop(0)
    
    # make a dictionary-entry per line
    for line in lines:
        values = line.rstrip().split()
        
        # make all values integers
        str_vals = values[1:]
        for i in range(len(str_vals)):
            str_vals[i] = int(str_vals[i])
            
        # add line to list
        ID_dict[int(values[0])] = str_vals
        
    return ID_dict


def normSig2_add(Days):
    """Preconditions:  Days is een lijst van dictionaries met per ID de 
                        bijbehorende waarden.
    Postconditions:    Normaliseert P2Sig (P2Sig * S1/S2) en voegt deze toe aan
                        elke dictionary."""
    for Day in Days:
        S1 = summation(Day,"P1Sig") # optellen van alle P1Sig waardes
        S2 = summation(Day,"P2Sig") # optellen van alle P2Sig waardes
        for ID in Day:
            Day[ID].append(int(Day[ID][Columns.index("P2Sig")])*S1/S2)


def summation(Day,column_head):
    """Preconditions:  Day is een dictionary met per ID de bijbehorende 
                        waarden en column head is de index van de kolom die
                        moet worden opgesomd.
    Postconditions:    Retourneert de som van alle waarden uit de kolom."""
    total = 0 # counter
    
    for ID in Day:
        total += int(Day[ID][Columns.index(column_head)])
        
    return total


def Log_add(Days):
    """Preconditions:  Days is een lijst van dictionaries met per ID de 
                        bijbehorende waarden.
    Postconditions:    Zet P1Sig, P2Sig en P2SigNorm in een logaritme en voegt 
                        deze toe aan elke dictionary in Days."""
    for Day in Days:
        for ID in Day:
            Day[ID].append(math.log10(Day[ID][Columns.index("P1Sig")]))
            Day[ID].append(math.log10(Day[ID][Columns.index("P2Sig")]))
            Day[ID].append(math.log10(Day[ID][Columns.index("P2SigNorm")]))

              
def Classes(Days,frequences=False):
    """Preconditions:  Days is een lijst van libraries, frequences is een 
                        boolean die het printen aan kan zetten, per default uit.
    Postconditions:    Roept de funtie make_class aan per dag en print daarmee
                        de frequenties per dag en klasse als frequences==True"""
    classes = "ABCD"
    
    for i in range(len(Days)):
        counts = make_class(Days[i])
        if frequences: # check whether we want the print output
            print("\nVoor dag",Day_numbers[i],"geldt: ")
            for j in range(len(counts)):
                print("\tFrequentie van",classes[j],"is",counts[j]) 


def make_class(Day,stb_threshold=25):
    """Preconditions:  Day is een dictionary met per ID de bijbehorende 
                        waarden. stb_threshold is de drempelwaarde voor de 
                        filtering van de classes.
    Postconditions:    Voegt aan elke waarde van Day een klasse (A, B, C of D)
                        toe op basis van de waarden van P1STB en P2STB. Voor
                        deze klassen geldt respectievelijk in vergelijking
                        met de drempelwaarde: Geen waarden lager, slechts P2 
                        lager, slechts P1 lager, beide lager.
                        Retourneert een lijst met de frequenties van alle
                        klassen"""
    counts = [0,0,0,0] # list of counters
    
    for ID in Day:
        if (Day[ID][Columns.index("P1STB")] >= stb_threshold):
            if (Day[ID][Columns.index("P2STB")] >= stb_threshold):
                Day[ID].append("A")
                counts[0]+=1
            else:
                Day[ID].append("B")
                counts[1]+=1
        else:
            if (Day[ID][Columns.index("P2STB")] >= stb_threshold):
                Day[ID].append("C")
                counts[2]+=1
            else:
                Day[ID].append("D")
                counts[3]+=1
                
    return counts   

  
def add_expression(Days):
    """Preconditions:  Days is een lijst van libraries
    Postconditions:    voegt de relatieve expressiewaarden toe aan elke Day uit
                        Days"""
    Columns.append("RelativeExpression")
    for Day in Days:
        for ID in Day:
            if (Day[ID][Columns.index("P1Sig")]/Day[ID][Columns.index("P2SigNorm")]) >= 1:
                r = (Day[ID][Columns.index("P1Sig")]/Day[ID][Columns.index("P2SigNorm")]) - 1
            else:
                r = (-Day[ID][Columns.index("P2SigNorm")]/Day[ID][Columns.index("P1Sig")]) + 1
            Day[ID].append(r)


def Daysdict(Days,r_filter=0.5):
    """Preconditions:  Days is een lijst van libraries, r_filter is de  
                        filtervoorwaarde voor de expressiewaarde r
    Postconditions:    Voegt alle r-waarden van Days samen in één dictionary 
                        en retourneert deze als eerste waarde. Houdt verder
                        bij van elk gen hoevaak deze elke belangrijke 
                        filtervoorwaarde overschrijdt, zodat hier later op 
                        gefilterd kan worden. Stopt deze waarden in een lijst
                        met respectievelijk de counters van: in class B of C,
                        in class D, spotgrootte buiten de range [40,160] en 
                        |r| kleiner dan r_filter."""
    # Dictionaries that will be filled
    New_Days = {}
    Filterinfo = {}
    
    for ID in Days[0]:
        New_Days[ID] = [] # the list of R-values that will become the value in the dictionary
        Filterinfo[ID] = [0,0,0,0]
        
        for i in range(len(Day_numbers)):
            New_Days[ID].append(Days[i][ID][-1]) # de laatste positie van de lijst kan worden gebruikt, omdat r altijd achteraan moet staan
            if -r_filter < Days[i][ID][-1] < r_filter: Filterinfo[ID][3] += 1 # telt hoe vaak |r| kleiner is dan de filterwaarde
            if not ((Days[i][ID][Columns.index("P1Cov")] >= 40) and (Days[i][ID][2] <= 160)): Filterinfo[ID][2] += 1 # telt hoe vaak de spotgrootte buiten de (arbitraire) waarden valt
            if (Days[i][ID][Columns.index("Class")] == "B") or (Days[i][ID][Columns.index("Class")] == "C"): Filterinfo[ID][0] += 1 # telt hoe vaak een gen in B of C zit
            if Days[i][ID][Columns.index("Class")] == "D": Filterinfo[ID][1] += 1 # telt hoe vaak een gen in D zit            
                
    return New_Days, Filterinfo

    
def filtering(rDict, Filterinfo):
    """Preconditions:  rDict is een library waarvan elke value een lijst is met
                        de r-waarde per dag. Filterinfo is een dictionary even 
                        groot als rDict, met per ID (key) een lijst (value)
                        met de hoeveelheid die elk gen de volgende
                        filtervoorwaarden respectievelijk overschrijdt: in 
                        class B of C, in class D, spotgrootte buiten de range 
                        [40,160] en |r| kleiner dan r_filter.
    Postconditions:    Retourneert een gefilterde versie van rDict op basis van
                        de filtervoorwaarden:
                            - D komt niet voor
                            - de r-waarde is minimaal éénmaal buiten de 
                              filtergrens
                            - de hoeveelheid keren dat C of B voorkomt, samen 
                              met hoe vaak de spotgrootte fout is mag niet meer 
                              zijn dan 2"""
    Filtered_rDict = {}
    
    for ID in rDict:
        if (Filterinfo[ID][0] + Filterinfo[ID][2] < 3) and (Filterinfo[ID][1] == 0) and (Filterinfo[ID][3] < 8): 
            Filtered_rDict[ID] = rDict[ID]
                
    return Filtered_rDict


def dict_to_txt(rDict,filename,parse):
    """Preconditions:  rDict is een dictionary van lijsten, filename is de naam
                        van het te schrijven bestand en parse is de ruimte die 
                        tussen elke waarde van rDict moet komen.
    Postconditions:    Schrijft een dictionary weg naar een tekstbestand met
                        één key-value paar per regel en een witregel aan het
                        einde"""
    outfile = open(filename,'w') # start writing in a file called filename
    
    for ID in rDict:
        outstring = str(ID)
        for r in rDict[ID]: outstring += parse + str('{:.3f}'.format(r))
        outfile.write(outstring+"\n")
        
    outfile.close()  


def plot_phase1(Days,rDict,r_filter=0.5):
    """Preconditions:  Days is een lijst van libraries
    Postconditions:    Roept alle mogelijkheden voor  
                        plot_Days aan."""   
    for log in [False,True]:
        for norm in [False,True]:
            plot_Days(Days,log,norm)
            
    plot_hist(rDict) # to plot the histograms
    line_plots(rDict,r_filter) # to plot the line plots
    

def plot_Days(Days,log=False,Norm=False):
    """Preconditions:  Days is een lijst van libraries, log en
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
    fig.suptitle("Visualisatie data"+add_on,size=24,weight='bold')
    for i in range(len(Days)):
        ax_day = ax[i%2,i//2] # in order to have the representation of the plots be plotted in 2 lines
        ax_day.set_title("Day "+str(Day_numbers[i]))
        df_day = pd.DataFrame.from_dict(Days[i]).transpose() # making a dataframe of rDict for more efficiency in plotting
        df_day.columns = Columns
        # assign color to right class for plotting
        df_day = df_day.replace({'Class': Colors})
        # plot every class per day; this ensures the legenda can be made
        for Class in Colors:
            # making a DataFrame for every class using a mask
            df_day_new = df_day[df_day['Class'] == Colors[Class]]
            # check whether we should use the logaritmic columns
            if log:
                # plotting one class
                df_day_new.plot(kind='scatter',x=Columns[7],y=Columns[8+norm_fac],c=Colors[Class],ax=ax_day,colorbar=False)
                # only plot P1 = P2 once, not four times
                if Class == 'A': ax_day.plot([1,5],[1,5],c="k")
            else: 
                # only plot P1 = P2 once, not four times
                if Class == 'A': ax_day.plot([0,50000],[0,50000],c="k")
                # plotting one class
                df_day_new.plot(kind='scatter',x=Columns[0],y=Columns[3+norm_fac],c=Colors[Class],ax=ax_day,colorbar=False)
            ax_day.legend(['P1=P2',"A","B","C","D"])

     
def plot_hist(rDict):
    """Preconditions:  rDict is een dictionary van lijsten,
    Postconditions:    Plot een histogram van de expressiewaarden voor elke 
                        dag, zodat er te zien is hoe vaak een expressiewaarde 
                        voorkomt."""
    fig, ax = plt.subplots(2,4,figsize=(20,10),sharex=True,sharey=True)
    fig.suptitle("Histogrammen expressiewaarden per dag",size=24,weight='bold')
    df = pd.DataFrame.from_dict(rDict).transpose() # making a dataframe of rDict for more efficiency in plotting

    for i in range(len(Day_numbers)):
        ax_day = ax[i%2,i//2]
        ax_day.set_title("Day "+str(Day_numbers[i]))
        ax_day.set_xlim(-5, 5)
        df[i].plot(kind='hist',ax=ax_day)


def line_plots(rDict,r_filter,movement=256):
    """Preconditions:  rDict is een dictionary van lijsten, r_filter is de 
                        filtergrens van |r| en movement is het aantal waarden
                        uit rDict die moeten worden overgeslagen voor het 
                        plotten, opdat andere plots gezien kunnen worden.
    Postconditions:    Plot ncols**2 lijndiagrammen van genen, met op elke x-as
                        de progressie van Days en op elke y-as de expressie-
                        waarde R. Plot tevens een lijn op y=r en y=-r om aan te 
                        geven waar de filtergrens zal liggen"""
    nr_of_plots = int(input("How many lineplots would you like? \n"))
    if nr_of_plots < 1: return # don't start plotting if the number of plots is 0 or negative
    
    ncols = math.ceil(nr_of_plots**0.5)
    # ncols = math.ceil(len(rDict)**0.5) # to plot all points (not recommended)
    nrows = math.ceil(nr_of_plots/ncols)
    fig, ax = plt.subplots(ncols,nrows,figsize=(ncols*5,ncols*5),squeeze=False,sharex=True,sharey=True)
    fig.suptitle(str(nr_of_plots)+" genexpressiewaarden door de tijd met filtergrenzen",size=ncols*10,weight='bold') 
    
    df_expr = pd.DataFrame.from_dict(rDict) # making a dataframe of rDict for more efficiency in plotting
    df_expr['Days'] = Day_numbers # adding the days in a column to make x-axis linear
    df_expr.set_index("Days",inplace=True)
    for col in range(nrows): # looping through 
        for row in range(ncols):
            if  col*ncols+row < nr_of_plots: # in order to only plot the amount asked
                ID = df_expr.columns[col*ncols+row+movement]
                df_expr[ID].plot(ax=ax[row,col])
                ax[row,col].set_title(str(ID))
                ax[row,col].axhline(y=r_filter,c='r') # plot upper line
                ax[row,col].axhline(y=-r_filter,c='r') # plot lower line
                ax[row,col].set_xticks(Day_numbers)
                ax[row,col].set_ylim(-2.5, 2.5)
    
    
def main(file_key_word):
    r_filter= 0.5                                           # to set the r_filter value
    outf_names = "Filtered_data.txt","Unfiltered_data.txt"
    Days = find_files(file_key_word)                               # make list of dictionaries
    normSig2_add(Days)                                     # normalise values
    Log_add(Days)                                          # calculate the LOG of values for plotting
    Classes(Days,frequences=True)                          # determine each values class
    add_expression(Days)                                   # calculate the r-value for expression per value
    rDict, Filterinfo = Daysdict(Days,r_filter)            # make one library for all r-values
    Filtered_rDict = filtering(rDict,Filterinfo)                     # filter the library based
    dict_to_txt(Filtered_rDict,outf_names[0],"  ")      # make a file of filtered data for phase 2
    dict_to_txt(rDict,outf_names[1],"  ")    # make a file of filtered data for phase 3
    if input("Would you like the plots of phase 1? [y]/n \n") == "y":    # ask whether the user would like to see the plots of phase 1 (boolean input)
        plot_phase1(Days,rDict,r_filter)                   # plot everything there is to plot in phase 1 based on unfiltered data
    
    return outf_names[0],outf_names[1]
