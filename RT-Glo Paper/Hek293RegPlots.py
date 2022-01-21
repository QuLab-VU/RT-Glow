#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 29 12:24:04 2021

@author: claytonwandishin
"""
#creates a reg plot from the HEK293 data
import os
import pandas as pd
from itertools import repeat
import datetime
import math
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn import preprocessing
from scipy import stats

hekdf = pd.read_csv('/Users/claytonwandishin/RT glow Model files/One_to_One_experiment/293dmsolumcell.csv')
hekdf=hekdf.drop(hekdf.index[22:79])
hekdf['cell.count'] = hekdf['cell.count'].astype(int)
hekdf['ch2.pos'] = hekdf['ch2.pos'].astype(int)
hekdf['lum']=hekdf['lum'].astype(int)
hekdf['time']=hekdf['time'].astype(float)
hekdf['drug1.conc']=hekdf['drug1.conc'].astype(int)
hekdf['Live_Dead']=hekdf['cell.count'] - hekdf['ch2.pos']
hekdf['log2_Live_Dead']=np.log2(hekdf['Live_Dead']+1)
hekdf['log_RLU']=np.log(hekdf['lum'])
cell_line_list = ['HEK-293FT']


#cell_line_list = completeDF['Cell_Line'].unique()
#cell_line_list = ['DMS53']
#drug_list=['DMSO']
for CL in cell_line_list:
    indivdata = hekdf.loc[hekdf['cell.line'] == CL]
    #indivdata = indivdata.loc[indivdata['Drug'] == 'hygromycin_b']
    #indivdata = indivdata.loc[indivdata['Well'] == 'G13']
    indivdrugdata = indivdata.loc[indivdata['drug1.conc'] == 0]
    #pd.to_numeric(indivdrugdata['TotHour_Lum'])
    indivdrugdata=indivdrugdata.loc[indivdrugdata['time'] < 100]
    #drops outliers
    indivdrugdata=indivdrugdata.drop([0,1,2,8])
    slope, intercept, r_value, pv, se = stats.linregress(indivdrugdata['log2_Live_Dead'],indivdrugdata['lum'])        
    sns.regplot(data=indivdrugdata, x=indivdrugdata['log2_Live_Dead'],    y=indivdrugdata['lum'], ci=95, truncate=False)
    #sns.regplot(data=indivdrugdata, x=indivdrugdata['TotHour'],    y=indivdrugdata['RLU'])
    sns.set(style='darkgrid')            
    #plt.xlim(0, 60)
    plt.ylabel(None)
    plt.xlabel(None)
    plt.yticks(size=15)
    plt.xticks(size=15)
    plt.legend(title=None, loc='upper left',labels=[r'$R^2:{0:.3f}$  CI:95%'.format(r_value)],prop={"size":16,"style":'oblique'})
    #plt.rcParams.update({'font.size': 25})
    #plt.title(CL)
    #plt.savefig('/Users/claytonwandishin/December 14 RT glow run/plots/'+CL+'NoDrug_LumCCReg100hrs.png', dpi=300,bbox_inches='tight')
    plt.show()
    print(r_value)
hekdf.to_csv('/Users/claytonwandishin/December 14 RT glow run/HEK293FTdf.csv')
