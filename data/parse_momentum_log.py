#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  4 15:12:06 2020

@author: claytonwandishin
"""
#parses the timestamps of reads from the momentum logs
import pandas as pd
from itertools import repeat

logDF = pd.read_csv('/Users/claytonwandishin/Downloads/momentum_log_CMW_11-03-2020.csv')

#truncated dataframe with just the columns with read, datetime, or barcode info
logDF_RDT = logDF.iloc[:,[5,11,13]]
logDF_RDT.columns = ['action','datetime','barcode']

#This dataframe now contains only the relevant information from the momentum log i.e. the rows where the action is read and the columns for datetime and barcode
logDF_RO = logDF_RDT.loc[logDF_RDT['action'] == 'Read']

#this creates a list to assign the "read number" to each barcode and timestamp pair by iterating through a range of 1-25 (our number of reads per plates) 6 times since we had 6 plates, this 

plates = [*range(1,26)]
read_number = []
for d in plates:   
    for i in repeat(None, 6):
       read_number.append(d)

logDF_RO['read_number'] = read_number

logDF_RO.to_csv('/Users/claytonwandishin/Dropbox/One to One RT Glow Experiment/20201103_LumFilesOnly/Parsed_MomentumLog.csv')