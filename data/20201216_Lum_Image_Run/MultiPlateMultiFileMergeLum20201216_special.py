#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 25 15:11:37 2020

@author: claytonwandishin
"""
#script for making a dataframe from the xls files
#needs to be run in the same folder that contains all the .xls file
import os
import pandas as pd
from itertools import repeat
import datetime
import math
import numpy as np

############################################
#MAKE SURE YOU ARE WORKING IN THE FOLDER WITH ALL THE FILES


##############################################
#list of cell lines in the experiment MAKE SURE THIS IS IN THE SAME ORDER AS THE PLATES WERE READ or not if you have the barcodes, but it'll make your life easier
cell_line_list = ['H1048','H841','DMS53','WM88']
#make sure these are correct
#number of plates
numplate = 4
#number of drug concnetrations in the series
drugser = 11
#number of reads per plate during the experiment
reads = 29
#this adds the drug info
extraplates=3
druglist1 = ['TAK-901','SCH1473759','YM-155','vemurafinib','TAK-901','SCH1473759','YM-155','vemurafinib','TAK-901','SCH1473759','YM-155','vemurafinib']
#multiply by number of drug concentrations
dlltotal = druglist1*drugser
dllextra = druglist1*drugser
druglist2 = ['AZD-1152','SNS-314','trametinib','hygromycin_b','AZD-1152','SNS-314','trametinib','hygromycin_b','AZD-1152','SNS-314','trametinib','hygromycin_b']
#multiply by number of drug concentrations
dll2 = druglist2*drugser
dlltotal.extend(dll2)
#multiply by number of plates
dlltotal = dlltotal*numplate
#multiply by number of reads
dlltotal = dlltotal*reads
dllextra.extend(dll2)
dllextra3 = dllextra*extraplates
dlltotal.extend(dllextra3)
dlltotal = pd.DataFrame(dlltotal)
dlltotal = dlltotal.rename(columns={0:'Drug'})


cll = cell_line_list
#multiply by number of reads
cll = cll*reads
cll.append('H1048')
cll.append('H841')
cll.append('DMS53')
cll = pd.DataFrame(cll)
#make sure range is from 1 to (numplate+1)
numberofplates = [*range(1,(numplate+1))]

#list of barcodes ideally ordered in same way as the cell lines or you need to create a dictionary, aka in the order they were read
barcode_list = ['V004832B','V004833B','V004834B','V004831B',]
bll = barcode_list
#multiply by number of reads
bll = bll*reads
bll.append('V004832B')
bll.append('V004833B')
bll.append('V004834B')
bll = pd.DataFrame(bll)

barcode_list = pd.DataFrame(barcode_list)
barcode_list = barcode_list.rename(columns={0:'barcode'})
cell_barcode_dict = cell_line_list

cbl = cell_barcode_dict

cbl = pd.DataFrame(cbl)
cbl = cbl.rename(columns={0:'cell_line'})
cbl['barcode'] = barcode_list['barcode']

#this walks through the file names to get the timestamp from the name
filelist = os.listdir()
filelist.sort()

timestamplist = []

#appends the timestamp slice of the filename to a list
for file in filelist:
    timestamplist.append((file[-16:-4]))
timestamplist = [int(x) for x in timestamplist]
timestamplist.sort()
timestamplist = [str(x) for x in timestamplist]

timestampdf = pd.DataFrame(timestamplist)
#converts the dataframe to datetimes
timestampdf = pd.to_datetime(timestampdf[0], format='%y%m%d%H%M%S')

#this creates a single dataframe with the timestamps associatd with each cell line
Totaldf = timestampdf
Totaldf = pd.DataFrame(Totaldf)
Totaldf = Totaldf.rename(columns={0:'DateTime'})
Totaldf['Cell_Line']= cll
Totaldf['Plate_Name']= bll


BigDF = pd.DataFrame()
#reads in the .xls file as a pandas dataframe

for m in filelist:
    # change this to use a variable 'file'
    df = pd.read_excel(m)
    
    # change this to use a variable 'file'
    filename = [m]
    timestamplist2 = []
    for j in filename:
        timestamplist2.append((j[-16:-4]))
    timestamplist2 = timestamplist2*264
    timestampdf2 = pd.DataFrame(timestamplist2)
    #converts the dataframe to datetimes
    timestampdf2 = pd.to_datetime(timestampdf2[0], format='%y%m%d%H%M%S') 
    timestampdf2 =pd.DataFrame(timestampdf2)  
    timestampdf2 = timestampdf2.rename(columns={0:'DateTime'})
    #this selects only the portion of the dataframe related to luminescence values
    
    ########################################
    #MAKE SURE THIS SLICE IS CORRECT BY CHECKING df VARIABLE
    dftrunc = df.iloc[33:45,3:25]
    #this melts that portion in to a single column df
    dfmelt = pd.melt(dftrunc)
    dfmelt = dfmelt.rename(columns={'value':'RLU'})
    dfmelt = pd.DataFrame(dfmelt)
    #make sure this matches the rows used in the experiment
    PlateRows = ['C','D','E','F','G','H','I','J','K','L','M','N']
    #creates a list in a range of column numbers and converts them to strings 
    PlateColumns = ["%02d" % x for x in range(2,24)]
    #a=["%02d" % x for x in range(24)]
    PlateColumns = [str(x) for x in PlateColumns] 
    
    #creates a well list compatible with the melted dataframe
    Wells=[]
    for j in PlateColumns:
        for k in PlateRows:
            Wells.append(k+j)
    Wells = pd.DataFrame(Wells)
    Wells = Wells.rename(columns={0:'Well'})
    LumAndWell = dfmelt['RLU']
    LumAndWell = pd.DataFrame(LumAndWell)
    LumAndWell['Well'] = Wells['Well']
    LumAndWell['DateTime']  = timestampdf2['DateTime']
    LumWellDateTime = LumAndWell
    BigDF = BigDF.append(LumWellDateTime)

#this last step merges the dataframe with all the data from the files with the dataframe that created a dictionary between timestamp and cell line to create one single dataframe with all of the data
BigDF = pd.merge(BigDF, Totaldf)
BigDF['Drug']=dlltotal['Drug']

#adds in a Drug_Conc column
#This is more complicated because the platemap starts at 0 then goes to the max conc. then does a gradient so I started with integers then looped
DrugCon = ['0','0.00000996','0.00000249','0.000000623','0.000000156','0.0000000389','0.00000000934','0.00000000233','0.000000000623','0.000000000156','0.0000000000778']

DrugConMult =[]
for d in DrugCon:   
    for i in repeat(None, 12):
        DrugConMult.append(d)

#This is because the concentration gradient is repeated twice across every plate since the drugs are side by side on a 384
#25 is the number of reads per plate or I guess the number of timestamps divided by the number of plates
#The addition of six is to account for the 3 extra plates or one missing plate rather in this specific dataset
DrugConMult = DrugConMult*((len(numberofplates)*2*reads)+6)
#converts all the elements in to floats because why not and then makes a dataframe
DrugConMult = [float(x) for x in DrugConMult]
DrugConMult = pd.DataFrame(DrugConMult)
DrugConMult = DrugConMult.rename(columns={0:'Drug_Conc'})

#adds in a drug units column
DrugUnits = ['M']
DrugUnits =  DrugUnits * len(BigDF['Drug'])
DrugUnits = pd.DataFrame(DrugUnits)
DrugUnits = DrugUnits.rename(columns={0:'Drug_Units'})
BigDF['Drug_Conc'] = DrugConMult['Drug_Conc']
BigDF['Drug_Units'] = DrugUnits['Drug_Units']

#make sure this tzero is the time that the drug was added in this experiment
tzero = datetime.datetime(2020,12,16,14,0,0)
BigDFdtdif = (BigDF['DateTime'] - tzero)
BigDFdtdif = (BigDFdtdif / np.timedelta64(1, 's'))/3600
BigDFdtdif = pd.DataFrame(BigDFdtdif)
#BigDFdtdif = BigDFdtdif.total_seconds()
BigDFdtdif = BigDFdtdif.rename(columns={'DateTime':'TotHour_Lum'})
BigDF['TotHour_Lum'] = BigDFdtdif['TotHour_Lum']
BigDF = BigDF.rename(columns={'DateTime':'Datetime_Lum'})
#BigDF.to_csv('/Users/claytonwandishin/December 14 RT glow run/CombinedLuminescenceTOTAL.csv')

#cbl.to_csv('/Users/claytonwandishin/December 14 RT glow run/Cell_Barcode_dict.csv')
