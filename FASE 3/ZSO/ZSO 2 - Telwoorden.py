# -*- coding: utf-8 -*-
"""
Created on Thu May  6 13:45:05 2021

@author: 20202876
"""


def TelWoorden(filename):
    
    inFile = open(filename)
    raw_data = inFile.readlines()
    inFile.close()

    data = [[line.strip()] for line in raw_data[1::3]]
    
    desc_freq = {}
    
    for description in data:
        if description not in desc_freq:
            desc_freq[description] = 1
        else:
            desc_freq[description] += 1
            
    desc_freq = {keys: values for keys, values in sorted(desc_freq.items(), key=lambda item: item[1])}
    
    return desc_freq
