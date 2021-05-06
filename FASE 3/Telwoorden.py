# -*- coding: utf-8 -*-
"""
Created on Thu May  6 00:47:02 2021

@author: 20202139
"""

from collections import Counter

def Telwoorden(bestandsnaam):
    infile = open(bestandsnaam)
    inlines = infile.read().splitlines()
    infile.close()
    
    for line in range(len(inlines)):
        print(inlines[line].replace("\\", " "))
   
    print(Counter(inlines))

Telwoorden('GenDescription2.txt')
            