"""
Course code: 8QA05
group number: 9
members: Iris Almekinders and Niels van Noort
"""


# importing needed libraries
import matplotlib.pyplot as plt
import pandas as pd
import math
from os import listdir
from pathlib import Path

# global variables
Columns = ['P2SigNorm', 'Log_P1Sig', 'Log_P2Sig', 'Log_P2SigNorm', 'Class']
Day_numbers = []

# definitions
def find_files(file_key_word):
    """Preconditions:   file_key_word is a string that occurs as a part of a 
                        file name. There is at least one file present in the 
                        same folder as this .py file whose name starts with
                        file_key_word, followed by a number, ending with .txt.
    Postconditions:     Returns list of textfile names based on Day_numbers."""

    # creating a list of filenames
    arr = listdir(Path(__file__).parent.absolute())
    filenames = [filename for filename in arr if file_key_word in filename]
    
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
    """Preconditions:   Filename is the name of a .txt file that can be found
                        in the same folder as this function.
    Postconditions:     Returns a dictionary containing a list of all 
                        expression values for each ID."""
                        
    # read the text file
    infile = open(filename)
    lines = infile.readlines()
    infile.close()
    ID_dict = {}
    if len(Columns) == 5: 
        for col_head in reversed(lines.pop(0).split()[1:]):
            Columns.insert(0,col_head) # adds the headings to global var Columns, only if that hasn't been done yet
    else: 
        lines.pop(0)
    
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
    """Preconditions:   Days is a list of dictionaries containing the expression
                        values for each ID.
    Postconditions:     Normalizes P2Sig (P2Sig * S1/S2) and adds this to every
                        dictionary."""
                        
    for Day in Days:
        S1 = summation(Day,"P1Sig") # counts all P1Sig values
        S2 = summation(Day,"P2Sig") # counts all P2Sig values
        for ID in Day:
            Day[ID].append(int(Day[ID][Columns.index("P2Sig")])*S1/S2)


def summation(Day,column_head):
    """Preconditions:   Day is a dictionary containing the expression values for
                        each ID and the column head is the index of the column
                        that has to be recited.
    Postconditions:     Returns the sum of all words of the column."""
    
    total = 0 # counter
    
    for ID in Day:
        total += int(Day[ID][Columns.index(column_head)])
        
    return total


def Log_add(Days):
    """Preconditions:   Days is a list of dictionaries containing the expression
                        values for each ID.
    Postconditions:     Puts P1Sig, P2Sig and P2SigNorm in a logarithm and adds
                        them to each dictionary in Days."""
                        
    for Day in Days:
        for ID in Day:
            Day[ID].append(math.log10(Day[ID][Columns.index("P1Sig")]))
            Day[ID].append(math.log10(Day[ID][Columns.index("P2Sig")]))
            Day[ID].append(math.log10(Day[ID][Columns.index("P2SigNorm")]))

              
def Classes(Days,frequences=False):
    """Preconditions:   Days is a list of libraries, frequences is a boolean that
                        can activate the printing, this is turned off by default.
    Postconditions:     Invokes the function make_class for every day and prints the
                        frequences for each day and class as frequences==True."""
        
    classes = "ABCD"
    
    for i in range(len(Days)):
        counts = make_class(Days[i])
        if frequences: # check whether we want the print output
            print("\nFor day",Day_numbers[i],"apply: ")
            for j in range(len(counts)):
                print("\tFrequency of",classes[j],"is",counts[j]) 


def make_class(Day,stb_threshold=25):
    """Preconditions:   Day is a dictionary containing the expression values for
                        each Id. stb_threshold is the treshold value for the
                        filtering of the classes.
    Postconditions:     Adds a class (A, B, C or D) to every value of Day based
                        on the values of P1STB and P2STB. For these values applies,
                        respectively in comparison with the treshold value,
                        No values lower, only P2 lower, only P1 lower, both lower.
                        Returns a list with the frequences of all classes."""
                        
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
    """Preconditions:   Days is a list of libraries.
    Postconditions:     Adds the relative expression values to each day in Days."""
    
    Columns.append("RelativeExpression")
    for Day in Days:
        for ID in Day:
            if (Day[ID][Columns.index("P1Sig")]/Day[ID][Columns.index("P2SigNorm")]) >= 1:
                r = (Day[ID][Columns.index("P1Sig")]/Day[ID][Columns.index("P2SigNorm")]) - 1
            else:
                r = (-Day[ID][Columns.index("P2SigNorm")]/Day[ID][Columns.index("P1Sig")]) + 1
            Day[ID].append(r)


def Daysdict(Days,r_filter=0.5):
    """Preconditions:   Days is a list of libraries, r_filter is the filter 
                        condition for expression value r.
    Postconditions:     Adds all r values of Days together in one dictionary
                        and returns this as the first value. It keeps track of
                        exceedances, to filter it later on. It puts these values
                        in a list with respectively the counters of:
                        in class B or C, in class D, spot size out of the range
                        [40, 160] and |r| smaller than r_filter."""
                        
    # Dictionaries that will be filled
    New_Days = {}
    Filterinfo = {}
    
    for ID in Days[0]:
        New_Days[ID] = [] # the list of R-values that will become the value in the dictionary
        Filterinfo[ID] = [0,0,0,0]
        
        for i in range(len(Day_numbers)):
            New_Days[ID].append(Days[i][ID][-1]) # the last position in the list can be used, because it's always the r-value
            if -r_filter < Days[i][ID][-1] < r_filter: Filterinfo[ID][3] += 1 # counts how often |r| is smaller than the filtervalue 
            if not ((Days[i][ID][Columns.index("P1Cov")] >= 40) and (Days[i][ID][2] <= 160)): Filterinfo[ID][2] += 1 # counts how often the spotsize is outside the (arbitrairy) values
            if (Days[i][ID][Columns.index("Class")] == "B") or (Days[i][ID][Columns.index("Class")] == "C"): Filterinfo[ID][0] += 1 # counts how often a gene is in class B or C
            if Days[i][ID][Columns.index("Class")] == "D": Filterinfo[ID][1] += 1 # counts how often a gene is in class D  
                
    return New_Days, Filterinfo

    
def filtering(rDict, Filterinfo):
    """Preconditions:   rDict is a library containing a list for each value, 
                        the list containing the r-value for each day. Filterinfo 
                        is a dictionary as big as rDict, containing a list 
                        (value) for each ID (key) with the amount that is 
                        exceeded for each gene: in class B or C, in class D, 
                        spot size outside the range [40, 160] and |r| smaller 
                        than r_filter.
    Postconditions:     Returns a filtered version of rDict based on the filter
                        conditions:
                            - D does not occur
                            - the r value is at least once outside the filter
                            limit
                            - the amount of times that C or B occur, together with
                            how often the spot size is wrong; this can't be more 
                            than 2"""
    Filtered_rDict = {}
    
    for ID in rDict:
        if (Filterinfo[ID][0] + Filterinfo[ID][2] < 3) and (Filterinfo[ID][1] == 0) and (Filterinfo[ID][3] < 8): 
            Filtered_rDict[ID] = rDict[ID]
                
    return Filtered_rDict


def dict_to_txt(rDict,filename,parse):
    """Preconditions:   rDict is a dictionary of lists, filename is the name of
                        the file that will be written, parse is the space 
                        between each value of rDict.
    Postconditions:     Makes a file from a dictionary with one key value pair 
                        for each line and a blank line at the end."""
                        
    outfile = open(filename,'w') # start writing in a file called filename
    
    for ID in rDict:
        outstring = str(ID)
        for r in rDict[ID]: outstring += parse + str('{:.3f}'.format(r))
        outfile.write(outstring+"\n")
        
    outfile.close()  


def plot_phase1(Days,rDict,r_filter=0.5):
    """Preconditions:  Days is a list of libraries.
    Postconditions:    Calls on all the possibilities for plot_Days.""" 
                        
    for log in [False,True]:
        for norm in [False,True]:
            plot_Days(Days,log,norm)
            
    plot_hist(rDict) # to plot the histograms
    line_plots(rDict,r_filter) # to plot the line plots
    

def plot_Days(Days,log=False,Norm=False):
    """Preconditions:   Days is a list from libraries, log and norm are booleans 
                        which indicate if the function must be plotted 
                        normalized or logarithmic.
    Postconditions:     Ensures that the data of every day would be in one 
                        figure.""" 
                        
    # make extra variables 
    add_on = '' # for the title 
    norm_fac = 0 # for normalisation column index
    Colors = {'A': 'red', 'B': 'blue', 'C': 'green', 'D': 'yellow'}
    
    # use the booleans for the title
    if log:
        add_on += ', logaritmic'
        if Norm:
            add_on += ' and normalized'
            norm_fac = 1  
    elif Norm: 
        add_on += ', normalized'
        norm_fac = 3
    
    # start plotting
    fig, ax = plt.subplots(2,4,figsize=(20,10),sharex=True,sharey=True)
    fig.suptitle("Visualisation data"+add_on,size=24,weight='bold')
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
    """Preconditions:   rDict is a dictionary of lists
    Postconditions:     plots a histogram from the expressionvalues for
                        each day, to visualize how often a particular
                        expression value occurs"""
                        
    fig, ax = plt.subplots(2,4,figsize=(20,10),sharex=True,sharey=True)
    fig.suptitle("Histograms of expressions per day",size=24,weight='bold')
    df = pd.DataFrame.from_dict(rDict).transpose() # making a dataframe of rDict for more efficiency in plotting

    for i in range(len(Day_numbers)):
        ax_day = ax[i%2,i//2]
        ax_day.set_title("Day "+str(Day_numbers[i]))
        ax_day.set_xlim(-5, 5)
        ax_day.set_ylabel("expression R")
        df[i].plot(kind='hist',ax=ax_day)


def line_plots(rDict,r_filter,movement=256):
    """Preconditions:   rDict is a dictionary of lists, r_filter is the filter
                        limit of |r|, movement is the amount of values
                        from rDict that need to be saved, so the other
                        plots can be seen
    Postconditions:     plots ncols**2 line charts of the genes. The x-axis
                        shows the progression of the days and the y-axis
                        the exression value R. It also plots the line y=r
                        and y=-r to indicate the boundaries"""
                        
    nr_of_plots = int(input("How many lineplots would you like? \n"))
    if nr_of_plots < 1: return # don't start plotting if the number of plots is 0 or negative
    
    ncols = math.ceil(nr_of_plots**0.5)
    # ncols = math.ceil(len(rDict)**0.5) # to plot all points (not recommended)
    nrows = math.ceil(nr_of_plots/ncols)
    fig, ax = plt.subplots(ncols,nrows,figsize=(ncols*5,ncols*5),squeeze=False,sharex=True,sharey=True)
    fig.suptitle(str(nr_of_plots)+" genexpression over time with filter boundaries",size=ncols*8,weight='bold') 
    
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
                ax[row,col].set_ylabel("Expression R")
    
    
def main(file_key_word):
    r_filter= 0.5                                                       # to set the r_filter value
    outf_names = "Filtered_data.txt","Unfiltered_data.txt"
    Days = find_files(file_key_word)                                    # make list of dictionaries
    normSig2_add(Days)                                                  # normalise values
    Log_add(Days)                                                       # calculate the LOG of values for plotting
    Classes(Days,frequences=True)                                       # determine each values class
    add_expression(Days)                                                # calculate the r-value for expression per value
    rDict, Filterinfo = Daysdict(Days,r_filter)                         # make one library for all r-values
    Filtered_rDict = filtering(rDict,Filterinfo)                        # filter the library based
    dict_to_txt(Filtered_rDict,outf_names[0],"  ")                      # make a file of filtered data for phase 2
    dict_to_txt(rDict,outf_names[1],"  ")                               # make a file of filtered data for phase 3
    if input("Would you like the plots of phase 1? [y]/n \n") == "y":   # ask whether the user would like to see the plots of phase 1 (boolean input)
        plot_phase1(Days,rDict,r_filter)                                # plot everything there is to plot in phase 1 based on unfiltered data
    
    return outf_names[0],outf_names[1]
