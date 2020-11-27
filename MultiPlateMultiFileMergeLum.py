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

#list of cell lines in the experiment MAKE SURE THIS IS IN THE SAME ORDER AS THE PLATES WERE READ
cell_line_list = ['293FT','WM88','WM1799','DMS53','H841','H1048']

#adds drug info to the csv
druglist1 = ['TAK-901','TAK-901','TAK-901','SCH1473759','SCH1473759','SCH1473759','Trametinib','Trametinib','Trametinib','SNS-314','SNS-314','SNS-314']
dlltotal = druglist1*11
druglist2 = ['YM-155','YM-155','YM-155','AZD-1152','AZD-1152','AZD-1152','Pimasertib','Pimasertib','Pimasertib','GSK-2879552','GSK-2879552','GSK-2879552']
dll2 = druglist2*11
dlltotal.extend(dll2)
dlltotal = dlltotal*6
dlltotal = dlltotal*25
dlltotal = pd.DataFrame(dlltotal)
dlltotal = dlltotal.rename(columns={0:'Drug'})


cll = cell_line_list
cll = cll*25
cll = pd.DataFrame(cll)
numberofplates = [*range(1,7)]

filelist = os.listdir()
filelist.sort()

timestamplist = []

#appends the timestamp slice of the filename to a list
for file in filelist:
    timestamplist.append((file[-16:-5]))
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



BigDF = pd.DataFrame()
#reads in the .xls file as a pandas dataframe

for m in filelist:
    # change this to use a variable 'file'
    df = pd.read_excel(m)
    
    # change this to use a variable 'file'
    filename = [m]
    timestamplist2 = []
    for j in filename:
        timestamplist2.append((j[-16:-5]))
    timestamplist2 = timestamplist2*264
    timestampdf2 = pd.DataFrame(timestamplist2)
    #converts the dataframe to datetimes
    timestampdf2 = pd.to_datetime(timestampdf2[0], format='%y%m%d%H%M%S') 
    timestampdf2 =pd.DataFrame(timestampdf2)  
    timestampdf2 = timestampdf2.rename(columns={0:'DateTime'})
    #this selects only the portion of the dataframe related to luminescence values
    dftrunc = df.iloc[53:65,3:25]
    #this melts that portion in to a single column df
    dfmelt = pd.melt(dftrunc)
    dfmelt = dfmelt.rename(columns={'value':'RLU'})
    dfmelt = pd.DataFrame(dfmelt)
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
BigDF.to_csv('/Users/claytonwandishin/Dropbox/One to One RT Glow Experiment/20201103_LumFilesOnly/CombinedLuminescence.csv')
