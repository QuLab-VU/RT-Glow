#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 22 14:53:29 2021

@author: claytonwandishin
"""
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
import scipy.optimize as opt
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy.stats import *

def errplot(x, y, yerr, **kwargs):
    ax = plt.gca()
    data = kwargs.pop("data")
    data.plot(x=x, y=y, yerr=yerr, kind="bar", ax=ax, **kwargs)
#for making a scatterplot from IC50 values


#This needs to be edited to pull in correct data
df1 = pd.read_csv('/Users/claytonwandishin/December 14 RT glow run/plots/PaperFigs/DipDimCompPlots/IC50scatterTidy.csv')
df1=df1.loc[df1['Drug']!='vemurafenib']
errordata1= pd.read_csv('/Users/claytonwandishin/December 14 RT glow run/plots/PaperFigs/DipDimCompPlots/ParamErrorDF.csv')
Ndata = pd.read_csv('/Users/claytonwandishin/December 14 RT glow run/plots/PaperFigs/DipDimCompPlots/Ndf.csv')
Ndata = Ndata[['Drug','Cell Line','Count Type','N']]
errordata2 = errordata1[['SD', 'Drug','Cell Line','Count Type','Parameter']]

errordata = errordata2.merge(Ndata)

errordata['SE']=errordata['SD']/np.sqrt(errordata['N'])
IC50err = errordata.loc[errordata['Parameter']=='EC50']
IC50err = IC50err.drop(columns='Parameter')
IC50err = IC50err.loc[IC50err['Drug']!='vemurafenib']
IC50err = IC50err.reset_index()
df1 = df1.merge(IC50err)
df1wil = df1.drop([df1.index[4],df1.index[16]])

sns.set(style='darkgrid')

df1 =df1.loc[df1['Cell Line']=='H1048']

#df1 = df1.loc[df1['Drug']!='TAK-901']
'''
dfdirect1=df1wil.loc[df1['Count Type']=='Direct']
dfdirect=dfdirect1['IC50']
dflum1=df1wil.loc[df1['Count Type']=='Lum']
dflum=dflum1['IC50']
wilcox=wilcoxon(dfdirect,dflum)
jj= wilcox.pvalue
mm =wilcox.statistic
print('pvalue = '+str(wilcox.pvalue))
'''
#plot1 = sns.catplot(data = df1, x='Drug', y='IC50',hue='Count Type',col='Cell Line', kind='strip')
#plot1.map(plt.errorbar, "Drug", "IC50", "SE")

ax = sns.pointplot('Drug', 'logIC50', hue='Count Type',data=df1, dodge=True, join=False,markers='_', ci=None, legend=None)


# Find the x,y coordinates for each point
x_coords = []
y_coords = []
for point_pair in ax.collections:
    for x, y in point_pair.get_offsets():
        x_coords.append(x)
        y_coords.append(y)

# Calculate the type of error to plot as the error bars
# Make sure the order is the same as the points were looped over
df1['95CI'] = df1['SE']*1.96

#logerror = np.log(df1['SE'])
#plt.errorbar(x_coords, y_coords, yerr=df1['95CI'], fmt='none',ecolor='k', elinewidth=1)




ax = sns.pointplot('Drug', 'logIC50', hue='Count Type',data=df1, dodge=True, join=False,markers='_', ci=None,legend=None)
plt.xticks(rotation = 45)
#plt.yscale('log')
#plt.ylim(-20,-3)
plt.ylabel('')
plt.xlabel('')
plt.legend('')
#plt.savefig('/Users/claytonwandishin/December 14 RT glow run/plots/PaperFigs/DipDimCompPlots/IC50comparisonH1048', dpi=300, bbox_inches='tight')
plt.show()
