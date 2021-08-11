#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 17 12:26:07 2021

@author: claytonwandishin

Modified by DRT 2021-04-27
"""
import os
import pandas as pd
from itertools import repeat
import datetime
import math
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

import os, sys
sys.path.insert(0, os.path.abspath('../'))


TOPDIR = [x for x in sys.path if 'RT-Glow' in x][0]
DATADIR = os.path.join(TOPDIR,'data')

#this reads in the re-ordered pyseg output csv file
cc = pd.read_csv(os.path.join(DATADIR,'cellcount20201214reorder.csv'))

cc.reset_index()

#these next lines create series from summing ever two rows of the data or just taking the value from every other row in the wellcond case, this will work as long as the imaging tiles == 2 if it is any other number then the value will need to be changed

dd = pd.DataFrame()
ccgb = cc['cell_count'].groupby(cc.index // 2).sum()
gfpgb = cc['ch2_pos'].groupby(cc.index // 2).sum()
wellcond = cc.iloc[::2, 6]
wellcond = wellcond.reset_index()

logDF = pd.read_csv(os.path.join(DATADIR,'MOMENTUM_LOGS_20201216.csv'))

#truncated dataframe with just the columns with read, datetime, or barcode info

logDF_RDT = logDF.iloc[:,[5,6,11,13]]
logDF_RDT.columns = ['action','read_type','datetime','barcode']

#This dataframe now contains only the relevant information from the momentum log i.e. the rows where the action is read and the columns for datetime and barcode
logDF_RO = logDF_RDT.loc[logDF_RDT['action'] == 'Read']
logDF_RO.reset_index()

ReadTime = logDF_RO.iloc[::2, 2]
ReadTime.reset_index()

#since there are 264 reads per individual plate, this creates a list going over each read time 264 times
readtime_plate = []
for d in ReadTime:   
    for i in repeat(None, 264):
       readtime_plate.append(d)
#since there are 264 reads per individual plate, this creates a list going over each barcode 264 times
cell_line_list = ['WM88','H1048','H841','DMS53']
cll_264 = []
for c in cell_line_list:   
    for i in repeat(None, 264):
       cll_264.append(c)
cll_264=cll_264*30

#since there are 264 reads per individual plate, this creates a list going over each barcode 264 times
barcode_list = ['V004831B','V004832B','V004833B','V004834B']    
bll_264 = []
for b in barcode_list:   
    for i in repeat(None, 264):
       bll_264.append(b)
bll_264=bll_264*30       

read = [*range(1,31)]
read_264_4plates =[]
for r in read:   
    for i in repeat(None, 1056):
       read_264_4plates.append(r)
#there are 264 wells per plate
ccDF = pd.DataFrame()
ccDF['Cell_Count'] = ccgb
ccDF['ch2_pos'] = gfpgb
ccDF['Live_Dead']= ccDF['Cell_Count']-ccDF['ch2_pos']
ccDF['Well']=wellcond['well']
ccDF['Cell_Line']=cll_264
ccDF['Plate_Name']=bll_264
ccDF['Read_Number']=read_264_4plates
ccDF['DateTime_Image']=readtime_plate

ccDF['DateTime_Image']=pd.to_datetime(ccDF['DateTime_Image'], format='%m/%d/%y %H:%M')
tzeroimage = datetime.datetime(2020,12,16,14,0,0)
ccDFdtdif = (ccDF['DateTime_Image'] - tzeroimage)
ccDFdtdif = (ccDFdtdif / np.timedelta64(1, 's'))/3600
ccDFdtdif = pd.DataFrame(ccDFdtdif)
ccDFdtdif =ccDFdtdif.rename(columns={'DateTime_Image':'TotHour_Image'})
ccDF['TotHour_Image']=ccDFdtdif['TotHour_Image']

output_filename = os.path.join(DATADIR,'CellCountDF_20201216.csv')
if !os.path.exists(output_filename):
    ccDF.to_csv(output_filename)
