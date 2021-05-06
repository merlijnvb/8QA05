# -*- coding: utf-8 -*-
"""
Created on Tue May  4 08:41:37 2021

@author: 20203520
"""


# def Bestand_inlezen(bestandsnaam):
#     infile = open(bestandsnaam)
#     inlines = infile.readlines()
#     infile.close()

# Bestand_inlezen('Voorbeeld_clusterresult.txt')
# Bestand_inlezen('GenDescription2.txt')

# def TelWoorden(bestandsnaam):
#     infile = open(bestandsnaam)
#     inlines = infile.read()
#     infile.close()
    
#     #print(inlines)
#     dic={}
#     inlines.lower()
#     for i in inlines.split()[1::3]:
#         if i in dic:
#             dic[i]= dic[i]+1
#         else:
#             dic[i]=1
#     print(dic)
     
def TelWoorden(bestandsnaam):
    infile = open(bestandsnaam)
    inlines = infile.readlines()
    infile.close()
    
    #print(inlines)
    dic={}
    for line in inlines[1::3]:
        line=line.lower()
        line=line.replace('(','')
        line=line.replace(')','')
        words=line.split()
        for i in words:
            if i in dic:
                dic[i]= dic[i]+1
            else:
                dic[i]=1
    print(dic)         
    print(len(dic.keys()))

TelWoorden('GenDescription2.txt')
    
    