#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 15 12:05:09 2021

@author: claytonwandishin
"""
import os
import pandas as pd
from itertools import repeat
import datetime
import math
import numpy as np
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
from sklearn import preprocessing
from scipy import stats
import scipy.optimize as opt



def ll4(x,b,c,d,e):
    '''This function is basically a copy of the LL.4 function from the R drc package with
     - b: hill slope
     - c: min response
     - d: max response
     - e: EC50'''
    return(c+(d-c)/(1+np.exp(b*(np.log(x)-np.log(e)))))
def ll3u(x,b,c,e):
    '''This function is basically a copy of the LL.3u function from the R drc package with
     - b: hill slope
     - c: min response
     - e: EC50'''
    #return(c+((1-c)/(1+(np.log(x/e)**b))))
    return(c+((1-c)/(1+np.exp(b*(np.log(x)-np.log(e))))))
    #plogIC50
    #return(c+((1-c)/(1+np.exp(b*(np.log(e)-np.log(x))))))
def pDose(x):
    '''This is just a helper function, to compute easily log transformed concentrations used in drug discovery'''
    return(np.log10(x))

def IC50(EC50,hill,minv):
    return(np.exp(np.log(EC50)+(1/hill)*np.log(1/(1-(2*minv)))))
#this script creates dual DIP/DIM figures from tidy data (separate df)

dipdf = pd.read_csv('/Users/claytonwandishin/December 14 RT glow run/HTS031_NCI_DIP_Viability_Relevant.csv')

dipdf['Count_Type']='Direct'

lumdf = pd.read_csv('/Users/claytonwandishin/December 14 RT glow run/20201216_Lum_CellCounts_TOTAL_DIM.csv')
lumdf = lumdf[['cell.line','drug1','drug1.conc','DIMr']]
lumdf['Count_Type']='Lum'
lumdf = lumdf.rename(columns={'DIMr':'dip_rate'})

lumzero = lumdf.loc[lumdf['drug1.conc'] == 0]
lumzerocll = lumzero['cell.line'].unique()
lumzerodll = lumzero['drug1'].unique()
lumzeroDIM = pd.DataFrame(columns=['cell.line','drug1','DIPzero'])
for c in lumzerocll:
    for d in lumzerodll:
        zerodimcalc = lumzero.loc[lumzero['cell.line'] == c]
        zerodimcalc = zerodimcalc.loc[zerodimcalc['drug1']== d]
        zeroDIMmean = zerodimcalc['dip_rate'].mean()
        lumzeroDIM=lumzeroDIM.append({"cell.line":c,"drug1":d,"DIPzero":zeroDIMmean}, ignore_index=True)

lumdf = lumdf.merge(lumzeroDIM)
lumdf['DIPrNorm']=lumdf['dip_rate']/lumdf['DIPzero']

dipzero = dipdf.loc[dipdf['drug1.conc'] == 0]
dipzerocll = dipzero['cell.line'].unique()
dipzerodll = dipzero['drug1'].unique()
dipzeroDIP = pd.DataFrame(columns=['cell.line','drug1','DIPzero'])
for c in dipzerocll:
    for d in dipzerodll:
        zerodipcalc = dipzero.loc[dipzero['cell.line'] == c]
        zerodipcalc = zerodipcalc.loc[zerodipcalc['drug1']== d]
        zeroDIPmean = zerodipcalc['dip_rate'].mean()
        dipzeroDIP=dipzeroDIP.append({"cell.line":c,"drug1":d,"DIPzero":zeroDIPmean}, ignore_index=True)

dipdf = dipdf.merge(dipzeroDIP)
dipdf['DIPrNorm']=dipdf['dip_rate']/dipdf['DIPzero'] 
cld = {'NCIH1048':'H1048','NCIH841':'H841'}    
dipdf = dipdf.replace({'cell.line': cld})        
lumdf = lumdf[lumdf['cell.line'] != 'DMS53']
colstokeep = list(lumdf.columns.values)
dipdf=dipdf[['cell.line','drug1','drug1.conc','dip_rate','Count_Type','DIPzero','DIPrNorm']]
lumdf = lumdf.drop_duplicates()


#lumdata needs to be fixed to only show the unique concentrations
dipdf = dipdf[dipdf['drug1'] != 'etoposide']
dipdf = dipdf[dipdf['drug1'] != 'AMG-900']
lumdf = lumdf[lumdf['drug1'] != 'hygromycin_b']


lumdrugs = lumdf['drug1'].unique()
#dipdrugs = dipdf['drug1'].unique()
#dipdrugs no barasertib
#dipdrugs = ['vemurafenib','trametinib','SNS-314','TAK-901','YM-155','SCH-1473759','barasertib']
dipdrugs=['barasertib']
dipdf=dipdf.append(lumdf)


#dipdf.to_csv('/Users/claytonwandishin/December 14 RT glow run/KScomparisonDF.csv')
#celllines = ['H841']
celllines = dipdf['cell.line'].unique()
errordf = pd.DataFrame()
dipdf = dipdf.rename(columns={'drug1.conc':'Drug_Con'})
dipdf['DClog']=pDose(dipdf['Drug_Con']+.0000000000001)
ll4paramlist=['hill coefficient','min response','max response','EC50']
ll3uparamlist=['hill coefficient','min response','EC50']
Ndf = pd.DataFrame(columns=['Cell Line','Count Type','Drug','N'])
CTlist=['Direct','Lum']

#this section creates a df with all of the N values for correct calculation of Standard Error
'''
for CL in celllines:
    for DD in dipdrugs:
        for ct in CTlist:
            workdf = dipdf.loc[dipdf['cell.line'] == CL]
            #this lines removes the zero values from the N calculation as the later fit does not use zeroes and therefore should not be included in the error calculation
            workdf = workdf.loc[workdf['Drug_Con'] != 0]
            workdf = workdf.loc[workdf['drug1'] == DD]
            workdf = workdf.loc[workdf['Count_Type'] == ct]
            n = len(workdf['cell.line'])
            Ndf=Ndf.append({'Cell Line':CL,'Drug':DD,'Count Type':ct,'N':n}, ignore_index=True)
  '''       
#Ndfprev = pd.read_csv('/Users/claytonwandishin/December 14 RT glow run/plots/PaperFigs/DipDimCompPlots/Ndf.csv')
#Ndf=Ndf.append(Ndfprev)
#Ndf.to_csv('/Users/claytonwandishin/December 14 RT glow run/plots/PaperFigs/DipDimCompPlots/Ndf.csv')
dipdf = dipdf.loc[dipdf['Drug_Con']!=0 ]

for CL in celllines:
    dipdfg = dipdf.loc[dipdf['cell.line'] == CL]
    dipdfg = dipdfg.loc[dipdfg['drug1']!= 'vemurafenib']
    dipdrugs=dipdfg['drug1'].unique()
    dipdrugs=['YM-155']
    for DD in dipdrugs:
        dipdfw = dipdf.loc[dipdf['cell.line'] == CL]
        dipdfw = dipdfw.loc[dipdfw['drug1'] == DD]
        #data fix for h841 and h1048 barasertib
        #dipdf = dipdf.loc[dipdf['DIPrNorm'] > -1.5]
        compoundData = dipdfw.groupby(['Count_Type'])
        fitData = []
        for name,group in compoundData:
            #ADDING BOUNDS TO THE FITTING FUNCTION IS CRUCIAL SO THAT LOG SCALED ITEMS DO NOT THROW A RUNTIME ERROR, WHICH GIVES MASSIVE SE
            fitCoefs, covMatrix = opt.curve_fit(ll3u, group.Drug_Con,group.DIPrNorm,bounds=([0,-5.,1e-12], [np.inf, 1.,np.inf]))
            resids = group.DIPrNorm-group.Drug_Con.apply(lambda x: ll3u(x,*fitCoefs))
            curFit = dict(zip(['b','c','e'],fitCoefs))
            curFit['Count_Type']=name
            curFit['residuals']=sum(resids**2)
            fitData.append(curFit)
            fitCompound = [ item['Count_Type'] for item in fitData]
            fitTable = pd.DataFrame(fitData).set_index('Count_Type')
            perr = np.sqrt(np.diag(covMatrix))
            errdf = pd.DataFrame(columns=['SD','SE','Cell Line','Drug','Count Type','Parameter'])
            cell = str(CL)
            durg = str(DD)
            errdf['SD']=perr
            #errdf['SE']=errdf['SD']/12
            errdf['Cell Line']= cell
            errdf['Drug']=durg
            errdf['Count Type']=name
            errdf['Parameter']=ll3uparamlist
            errordf = errordf.append(errdf)
            if name == 'Direct':
                IC50valDirect =       IC50(fitTable['e'],fitTable['b'],fitTable['c']).copy()
                IC50valDirect=IC50valDirect.reset_index()
                IC50valDirect2=IC50valDirect.at[0,0]
                IC50direct = r'IC50 : {:0.3e}'.format(IC50valDirect2)+' M'
            if name == 'Lum':
                IC50valLum =       IC50(fitTable['e'],fitTable['b'],fitTable['c']).copy()
                IC50valLum=IC50valLum.reset_index()
                IC50valLum2=IC50valLum.at[1,0]
                IC50lum =r'IC50 : {:0.3e}'.format(IC50valLum2)+' M'
            refDose = -np.linspace(min(dipdfw['DClog'])*0.8,max(dipdfw['DClog'])*3.9,256)
            refDose = (10**-refDose)*1e6
            #sns.lmplot(x='DClog',y='DIPrNorm',data=dipdfw,hue='Count_Type',fit_reg=False)._legend.set_title('Count Type')
            sns.lmplot(x='DClog',y='DIPrNorm',data=dipdfw,hue='Count_Type',fit_reg=False, legend=False)
            for fit in fitData:
                plt.plot([pDose(i) for i in refDose],[ll3u(i,*[fit[i] for i in ['b','c','e']]) for i in refDose])
        
       
        
                        
        sns.set(style='darkgrid')     
        '''
        plt.ylabel('Normalized DIP or Luminescence Rate')
        plt.xlabel('log Drug Concentration (M)')
        plt.title('CL+" treated with "+DD')
        '''
        plt.ylabel('')
        plt.xlabel('')
        plt.title('')
        plt.yticks(fontsize=20)
        plt.xticks(fontsize=20)
        #plt.legend(title=None, loc='upper right', labels=[IC50direct,IC50lum])
        plt.savefig('/Users/claytonwandishin/December 14 RT glow run/plots/PaperFigs/DipDimCompPlots/'+CL+" treated with "+DD, dpi=500, bbox_inches='tight')
        #plt.savefig('/Users/claytonwandishin/December 14 RT glow run/plots/PaperFigs/DipDimCompPlots/'+CL+" treated with "+DD, dpi=500)
        plt.show()
#edfprev =pd.DataFrame()
#edfprev = pd.read_csv('/Users/claytonwandishin/December 14 RT glow run/plots/PaperFigs/DipDimCompPlots/ParamErrorDF.csv')
#edfprev = edfprev.append(errordf)
#edfprev.to_csv('/Users/claytonwandishin/December 14 RT glow run/plots/PaperFigs/DipDimCompPlots/ParamErrorDF.csv')
