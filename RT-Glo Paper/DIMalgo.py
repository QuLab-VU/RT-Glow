#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 31 14:06:51 2021

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

#this script is for determing the slope of the luminescence or direct count data when accounting for minimum luminescence being reached, since if you don't it slopes up after the first total potency concentration
 
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
def pDose(x):
    '''This is just a helper function, to compute easily log transformed concentrations used in drug discovery'''
    return(np.log10(x))

def IC50(EC50,hill,maxv,minv):
    return(np.exp(np.log(EC50)+(1/hill)*np.log(maxv/(maxv-(2*minv)))))
debugdf = pd.read_csv('/Users/claytonwandishin/December 14 RT glow run/20201216_Lum_CellCounts_TOTAL.csv')
data = pd.read_csv('/Users/claytonwandishin/December 14 RT glow run/050819_RTGlowDF.csv')
#data =data.rename(columns={'cell.line':'Cell_Line','well':'Well', 'drug1':'Drug','drug1.conc':'Drug_Conc'})

#df.rename(columns={"A": "a", "B": "c"})
data = data.rename(columns={'TotHour':'TotHour_Lum', 'well':'Well', 'drug1':'Drug','cell.line':'Cell_Line', 'drug1.conc':'Drug_Conc'})
debugdf = debugdf.rename(columns={'TotHour':'TotHour_Lum', 'well':'Well', 'drug1':'Drug','cell.line':'Cell_Line', 'drug1.conc':'Drug_Conc'})

drugLL = data['Drug'].unique()


#####used during debug



#data = data[data.Well != 'J12']
#data = data[data.Well != 'K12']
#data = data[data.Well != 'K10']
#data = data[data.Well != 'K11']











##### THIS IS WHERE THE LUM RATES ARE STORED
dimrate = pd.DataFrame(columns=['Well','DIMr','Cell_Line','Drug','R2_DIM_Slice','DIM_Slice'])

IC50df = pd.DataFrame(columns=['Cell_Line','Drug','IC50_Lum'])
cline = ['H1930']
drug = ['SNS-314']
druglistref = list(data.Drug.unique())
#conc = [0, 9.96e-06,2.49e-06,6.23e-07]
#conc = [0,9.96e-06]
rrLIST = []
conc = list(data.Drug_Conc.unique())
conc.sort()
#cldata = data.loc[data['Cell_Line'] == 'CORL279']
#cldrugdata = cldata.loc[cldata['Drug'] == 'barasertib']
#for d in conc:
    #sns.lineplot(data=cldrugdata.loc[cldrugdata['Drug_Conc'] == d], x='TotHour_Lum', y='RLU')
    #plt.show()
dcrr = []
dcrrslice =[]
dcrrlabel =[]
dcrrwell=[]
for CL in cline:
    for dd in drug:
        for cc in conc:
            cldata = data.loc[data['Cell_Line'] == CL]
            cldrugdata = cldata.loc[cldata['Drug'] == dd]
            cldrugdata11c = cldrugdata.loc[cldrugdata['Drug_Conc'] == cc]
            welly = list(cldrugdata11c.Well.unique())
            for w in welly:
                cldrugdata1c = cldrugdata11c.loc[cldrugdata11c['Well'] == w]
                if cc == 0.0:
                    maxlum0 = cldrugdata1c['RLU'].max()
                    #maxlum0 = 120
                    #the addition of iloc at the end specifies it as a single value instead of what would be a series object but make sure this works otherwise you need to manually put in the maxtime to generate slopes from
                    #maxtime0 = cldrugdata1c.query('RLU == '+str(maxlum0))['TotHour_Lum'].iloc[0]
                    #maxtime CORL279 AMG900
                    #maxtime0=100
                    #maxtime for H1930 SCH1473759
                    maxtime0 = 93
                    #maxtime0 for H841 AZD1152 (barasertib)
                    #maxtime0= 117.82805555555555
                    #maxtime0 for CORL279 SCH1473759
                    #maxtime0 = 47.5536
                    #H1048 adjustment for TAK901
                    #maxtime0 = 70
                    #H1048 Barasertib adjustment
                    #maxtime0=120
                    #maxtime0=85

                    cldrugdata1c = cldrugdata1c.loc[cldrugdata1c['TotHour_Lum'] <= maxtime0]
                    cldrugdata1c=cldrugdata1c.loc[cldrugdata1c['Well'] == w]
                    '''
                    slope, intercept, r_value, pv, se = stats.linregress(cldrugdata1c['TotHour_Lum'],cldrugdata1c['RLU'])
                    rr = r_value**2
                    rrLIST.append(rr)
                    dimrate = dimrate.append({'Well':w,'DIMr':slope}, ignore_index=True)
                    '''
                    sns.set(style='darkgrid') 
                    sns.regplot(data=cldrugdata1c, x=cldrugdata1c['TotHour_Lum'],    y=cldrugdata1c['RLU'], ci=95, fit_reg=False)
 #############################

                    maxtime=maxtime0
                    cldrugdata3c = cldrugdata1c.loc[cldrugdata1c['TotHour_Lum'] <= maxtime]
                        #here a loop needs to be written to take the unique time points that are left and put them in a list, then from that list take increasing slices as ilocs from the end stored as a new dataframe until the R2 value reaches >.90, then store that time as an object mintime and slice the original dataframe from that mintime >=
                    timepoints = cldrugdata3c['TotHour_Lum'].unique()
                    numtp = len(timepoints)-3
                    #something gets messed up here sometimes but what needs to be accomlished is not allowing the algo to take a slope of a 3 point line
                    intloclisttp = [*range(0,numtp)]
                    rrlist2=[]
                    ltp = len(timepoints)
                    for t in intloclisttp:
                        slic=timepoints[t]
                        slicesize=ltp-t
                        mintime = slic
                        #print(str(mintime))
                        cldrugdata2c = cldrugdata3c.loc[cldrugdata3c['TotHour_Lum'] >= mintime]
                        slope, intercept, r_value, pv, se = stats.linregress(cldrugdata2c['TotHour_Lum'],cldrugdata2c['RLU'])
                        rr=r_value**2
                        dcrr.append(rr)
                        dcrrslice.append(slicesize)
                        dcrrlabel.append(cc)
                        dcrrwell.append(w)
                        rrlist2.append(rr)
                        #print(str(rr))
                    #CONVERT THIS TO A SERIES
                    rrlist3=pd.Series(rrlist2)
                    rrmax = rrlist3[rrlist3 == rrlist3.max()].index[0]
                    mintime=timepoints[rrmax]
                    cldrugdata1c = cldrugdata1c.loc[cldrugdata1c['TotHour_Lum'] >= mintime]
                    SlicePoints = list(cldrugdata1c['TotHour_Lum'])
                    slope, intercept, r_value, pv, se = stats.linregress(cldrugdata1c['TotHour_Lum'],cldrugdata1c['RLU'])
                    rr=r_value**2
                    dimrate = dimrate.append({'Well':w,'DIMr':slope,'Cell_Line':CL,'Drug':dd,'R2_DIM_Slice':rr,'DIM_Slice':SlicePoints}, ignore_index=True)
                    print('Maximum R2 for '+CL+' treated with '+dd+' at '+str(cc)+'in well '+w+'is '+str(rr))
                    sns.regplot(data=cldrugdata1c, x=cldrugdata1c['TotHour_Lum'],    y=cldrugdata1c['RLU'], ci=95)
                    #sns.regplot(data=indivdrugdata, x=indivdrugdata['TotHour'],    y=indivdrugdata['RLU'])
                    sns.set(style='darkgrid')            
                    #plt.xlim(0, 100)
                    plt.ylabel('')
                    plt.xlabel('')
                    plt.xlim(0,maxtime)
                    plt.legend(title=None, loc='upper left', labels=[r'$R^2:{0:.3f}$  CI:95%'.format(r_value**2)])
                    #plt.rcParams.update({'font.size': 25})
                    plt.title(CL+' treated with '+str(cc)+'M '+dd+' in well '+w)
                    plt.savefig('/Users/claytonwandishin/December 14 RT glow run/plots/PaperFigs/DRC_Lum/H1048/AlgoV2/'+CL+'_'+dd+'/IndivRegLines/'+w+'at'+str(cc)+'Reg'+'.png', dpi=300, bbox_inches='tight')
                    plt.show()










  ########                 
                    if rr <.90:
                        print('Regression model for '+CL+' treated with '+dd+' at '+str(cc)+'M is <.95 MODEL DOES NOT FIT R value ='+str(r_value**2))
                        sns.regplot(data=cldrugdata1c, x=cldrugdata1c['TotHour_Lum'],    y=cldrugdata1c['RLU'], ci=95)
                        #sns.regplot(data=indivdrugdata, x=indivdrugdata['TotHour'],    y=indivdrugdata['RLU'])
                        sns.set(style='darkgrid')            
    #plt.xlim(0, 100)
                        plt.ylabel('')
                        plt.xlabel('')
                        plt.xlim(0,maxtime0)
                        plt.legend(title=None, loc='upper left', labels=[r'$R^2:{0:.3f}$  CI:95%'.format(r_value**2)])
                        #plt.rcParams.update({'font.size': 25})
                        plt.title(CL+' treated with '+str(cc)+'M '+dd+' in '+w)
                        #plt.savefig('/Users/claytonwandishin/December 14 RT glow run/plots/PaperFigs/DRC_Lum/H1048/AlgoV2/H1048_Barasertib/'+w+'at'+str(cc)+'.png', dpi=300, bbox_inches='tight')
                        plt.show()
                    elif rr>.90:
                        print('Regression model for '+CL+' treated with '+dd+' at '+str(cc)+' is >.95 PASS')
                        sns.regplot(data=cldrugdata1c, x=cldrugdata1c['TotHour_Lum'],    y=cldrugdata1c['RLU'], ci=95)
                        #sns.regplot(data=indivdrugdata, x=indivdrugdata['TotHour'],    y=indivdrugdata['RLU'])
                        sns.set(style='darkgrid')            
    #plt.xlim(0, 100)
                        plt.ylabel('')
                        plt.xlabel('')
                        plt.xlim(0,maxtime0)
                        plt.legend(title=None, loc='upper left', labels=[r'$R^2:{0:.3f}$  CI:95%'.format(r_value**2)])
                        #plt.rcParams.update({'font.size': 25})
                        plt.title(CL+' treated with '+str(cc)+'M '+dd+' in '+w)
                        #plt.savefig('/Users/claytonwandishin/December 14 RT glow run/plots/PaperFigs/DRC_Lum/H1048/AlgoV2/H1048_Barasertib/'+w+'at'+str(cc)+'.png', dpi=300, bbox_inches='tight')
                        plt.show()
                        
                elif cc != 0:
                    #this minimun value may be cell line and or drug depndent and this code should be edited to reflect that
                    minlum=1000
                    cldrugdata1c =      cldrugdata1c.loc[cldrugdata1c['TotHour_Lum'] <= maxtime0]
                    cldrugdata1c=cldrugdata1c.loc[cldrugdata1c['Well'] == w]
                    mincheckRLUcol = cldrugdata1c['RLU'].min()
                    maxldrug = cldrugdata1c['RLU'].max()
                    #max time needs to first check if the minimum value is reached for that drug concentration, then if not it needs to use the maxtime value from the zero concentration and use that as the maxtime
                    #This plotting function shows ALL of the points without a regression which is necessary for pretty figures, by adding to the regplot below we get the orange regression line
                    sns.regplot(data=cldrugdata1c, x=cldrugdata1c['TotHour_Lum'],    y=cldrugdata1c['RLU'], ci=95, fit_reg=False)
                    #This should be adjusted to loop through like the else function below
                    if mincheckRLUcol <= minlum:               
                        maxtime = cldrugdata1c.query('RLU <= '+str(minlum))['TotHour_Lum'].iloc[0]
                        cldrugdata1c = cldrugdata1c.loc[cldrugdata1c['TotHour_Lum'] <= maxtime]
                        slope, intercept, r_value, pv, se = stats.linregress(cldrugdata1c['TotHour_Lum'],cldrugdata1c['RLU'])
                        rr=r_value**2
                        SlicePoints = list(cldrugdata1c['TotHour_Lum'])
                        dimrate = dimrate.append({'Well':w,'DIMr':slope,'Cell_Line':CL,'Drug':dd,'R2_DIM_Slice':rr,'DIM_Slice':SlicePoints}, ignore_index=True)
                        sns.regplot(data=cldrugdata1c, x=cldrugdata1c['TotHour_Lum'],    y=cldrugdata1c['RLU'], ci=95, fit_reg=False)
                        #sns.regplot(data=indivdrugdata, x=indivdrugdata['TotHour'],    y=indivdrugdata['RLU'])
                        sns.set(style='darkgrid')            
                        #plt.xlim(0, 100)
                        plt.ylabel('')
                        plt.xlabel('')
                        plt.xlim(0,maxtime)
                        plt.legend(title=None, loc='upper left', labels=[r'$R^2:{0:.3f}$  CI:95%'.format(r_value**2)])
                        #plt.rcParams.update({'font.size': 25})
                        plt.title(CL+' treated with '+str(cc)+'M '+dd+' in well '+w)
                        plt.savefig('/Users/claytonwandishin/December 14 RT glow run/plots/PaperFigs/DRC_Lum/H1048/AlgoV2/'+CL+'_'+dd+'/IndivRegLines/'+w+'at'+str(cc)+'.png', dpi=300, bbox_inches='tight')
                        plt.show()
                    else:
                        maxtime = maxtime0
                        mintime = 0
                        cldrugdata3c = cldrugdata1c.loc[cldrugdata1c['TotHour_Lum'] <= maxtime]
                        #here a loop needs to be written to take the unique time points that are left and put them in a list, then from that list take increasing slices as ilocs from the end stored as a new dataframe until the R2 value reaches >.90, then store that time as an object mintime and slice the original dataframe from that mintime >=
                        timepoints = cldrugdata3c['TotHour_Lum'].unique()
                        numtp = len(timepoints)-3
                        #something gets messed up here sometimes but what needs to be accomlished is not allowing the algo to take a slope of a 2 point line
                        intloclisttp = [*range(0,numtp)]
                        rrlist2=[]
                        ltp = len(timepoints)
                       
                        for t in intloclisttp:
                            slicesize=ltp-t
                            slic=timepoints[t]
                            mintime = slic
                            #print(str(mintime))
                            cldrugdata2c = cldrugdata3c.loc[cldrugdata3c['TotHour_Lum'] >= mintime]
                            slope, intercept, r_value, pv, se = stats.linregress(cldrugdata2c['TotHour_Lum'],cldrugdata2c['RLU'])
                            rr=r_value**2
                            dcrr.append(rr)
                            dcrrslice.append(slicesize)
                            dcrrlabel.append(cc)
                            dcrrwell.append(w)
                            rrlist2.append(rr)
                            #print(str(rr))
                        #CONVERT THIS TO A SERIES
                        rrlist3=pd.Series(rrlist2)
                        rrmax = rrlist3[rrlist3 == rrlist3.max()].index[0]
                        mintime=timepoints[rrmax]
                        cldrugdata1c = cldrugdata1c.loc[cldrugdata1c['TotHour_Lum'] >= mintime]
                        SlicePoints = list(cldrugdata1c['TotHour_Lum'])
                        slope, intercept, r_value, pv, se = stats.linregress(cldrugdata1c['TotHour_Lum'],cldrugdata1c['RLU'])
                        rr=r_value**2
                        dimrate = dimrate.append({'Well':w,'DIMr':slope,'Cell_Line':CL,'Drug':dd,'R2_DIM_Slice':rr,'DIM_Slice':SlicePoints}, ignore_index=True)
                        print('Maximum R2 for '+CL+' treated with '+dd+' at '+str(cc)+'in well '+w+'is '+str(rr))
                        sns.regplot(data=cldrugdata1c, x=cldrugdata1c['TotHour_Lum'],    y=cldrugdata1c['RLU'], ci=95)
                        #sns.regplot(data=indivdrugdata, x=indivdrugdata['TotHour'],    y=indivdrugdata['RLU'])
                        sns.set(style='darkgrid')            
                        #plt.xlim(0, 100)
                        plt.ylabel('')
                        plt.xlabel('')
                        plt.xlim(0,maxtime)
                        plt.legend(title=None, loc='upper left', labels=[r'$R^2:{0:.3f}$  CI:95%'.format(r_value**2)])
                        #plt.rcParams.update({'font.size': 32})
                        plt.title(CL+' treated with '+str(cc)+'M '+dd+' in well '+w)
                        plt.savefig('/Users/claytonwandishin/December 14 RT glow run/plots/PaperFigs/DRC_Lum/H1048/AlgoV2/'+CL+'_'+dd+'/IndivRegLines/'+w+'at'+str(cc)+'.png', dpi=300, bbox_inches='tight')
                        plt.show()
#dimrate = dimrate.fillna()
full_df = cldrugdata.copy()
full_df = full_df.merge(dimrate)
dimratezero = full_df.loc[full_df['Drug_Conc'] == 0]
dimratezeromean = pd.Series(dimratezero['DIMr'].unique()).mean()
dzmean = dimratezeromean
full_df['DIMrNorm'] = full_df['DIMr']/dzmean
#double check to make sure that plotting a datapoint for zero is the best way to do this if the data has already been normalized
#something is weird with the scaling
full_df['DClog']= np.log10(full_df['Drug_Conc']+.000000000001)
#full_df['DIMrNormal'] =
#sns.lmplot(x='DClog',y='DIMrNorm',data=full_df,fit_reg=False)
#removes erroneous points for the H1048 YM-155 data
#full_df = full_df[full_df.Drug_Conc != 9.340000000000001e-09]
#full_df = full_df[full_df.Drug_Conc != 3.89e-08]
#removes errod in H841 TAK901 dataset
#full_df = full_df[full_df.Drug_Conc != 9.96e-06]
#removes error in the H841trametinib dataset
#full_df = full_df[full_df.Drug_Conc != 7.779999999999999e-11]
#removes error in the H841 SNS314 dataset
#full_df = full_df[full_df.Drug_Conc != 7.779999999999999e-11]
#this well error could also be fixed by forcing the slope
#full_df = full_df[full_df.Well != 'D22']
#this well error could also be fixed by forcing slope H526 with barasertib
#full_df = full_df[full_df.Well != 'H11']
#this is well errors for DMS114 with barasertib
#full_df = full_df[full_df.Well != 'F11']
#full_df = full_df[full_df.Well != 'G11']
#full_df = full_df[full_df.Well != 'H11']
#full_df = full_df[full_df.Well != 'F07']
#full_df = full_df[full_df.Well != 'G07']
#full_df = full_df[full_df.Well != 'H07']
#DMS114 with AMG900
#full_df = full_df[full_df.Well != 'C12']
#fixes well error for DMS114TAK901
#full_df = full_df[full_df.Well != 'D20']
#fixes well error for H1048 with barasertib
#full_df = full_df[full_df.Well != 'K23']
#full_df = full_df[full_df.Well != 'K16']
#fixes H1048 YM-155 errors
#full_df = full_df[full_df.Well != 'E08']
#full_df = full_df[full_df.Well != 'I08']
#full_df = full_df[full_df.Well != 'M08']
#fixes h1048 trametinib
#full_df = full_df[full_df.Well != 'I22']
#fixes h1048 SNS314
#full_df = full_df[full_df.Well != 'L13']
#full_df = full_df[full_df.Well != 'D14']
#fixes H841 barasertib
#full_df = full_df[full_df.Well != 'K13']
#fixes H841 hygromycin
#full_df = full_df[full_df.Well != 'N23']
#fixes H841 SNS314
#full_df = full_df[full_df.Well != 'H22']
#before applying this, the entire dataframe needs to be truncated so that the residuals make more sense
#DMS1114 TAK901
#full_df = full_df[full_df.Well != 'E20']
#full_df = full_df[full_df.Well != 'E19']
#DMS114 YM155
#full_df = full_df[full_df.DIMrNorm < 1.24]
#full_df = full_df[full_df.Well != 'F18']
#full_df = full_df[full_df.Well != 'G18']
#full_df = full_df[full_df.Well != 'H18']
#CORL279 SCH1473759 fix
#full_df = full_df[full_df.Well != 'K10']
#full_df = full_df[full_df.Well != 'K11']
#full_df = full_df[full_df.Well != 'K12']
#H1930 SCH1473759
#full_df = full_df[full_df.DIMrNorm > -1.0]
#H1930 etoposide
'''
full_df = full_df[full_df.Well != 'J22']
full_df = full_df[full_df.Well != 'K22']
full_df = full_df[full_df.Well != 'J21']
full_df = full_df[full_df.Well != 'K21']
full_df = full_df[full_df.Well != 'J20']
full_df = full_df[full_df.Well != 'K20']
full_df = full_df[full_df.Well != 'I19']
full_df = full_df[full_df.Well != 'I18']
full_df = full_df[full_df.Well != 'K18']
full_df = full_df[full_df.Well != 'K14']
full_df = full_df[full_df.Well != 'K15']
full_df = full_df[full_df.Well != 'I15']
full_df = full_df[full_df.Well != 'K16']
full_df = full_df[full_df.Well != 'I16']
full_df = full_df[full_df.Well != 'J17']
'''
#H1930 SNS-314 data does not look good
'''
full_df = full_df[full_df.Well != 'L12']
full_df = full_df[full_df.Well != 'M12']
full_df = full_df[full_df.Well != 'N12']
full_df = full_df[full_df.Well != 'L11']
full_df = full_df[full_df.Well != 'M11']
full_df = full_df[full_df.Well != 'N11']
full_df = full_df[full_df.Well != 'M10']
full_df = full_df[full_df.Well != 'N10']
full_df = full_df[full_df.Well != 'N07']
full_df = full_df[full_df.Well != 'N03']
'''
compoundData = full_df.groupby(['Drug'])
fitData = []
for name,group in compoundData:
    fitCoefs, covMatrix = opt.curve_fit(ll4, group.Drug_Conc,         group.DIMrNorm)
    resids = group.DIMrNorm-group.Drug_Conc.apply(lambda x: ll4(x,*fitCoefs))
    curFit = dict(zip(['b','c','d','e'],fitCoefs))
    curFit['Drug']=name
    curFit['residuals']=sum(resids**2)
    fitData.append(curFit)
    fitCompound = [ item['Drug'] for item in fitData]
    fitTable = pd.DataFrame(fitData).set_index('Drug')
    
    #Export this fit table into a physical file
    
    
    
    
    IC50valLum =       IC50(fitTable['e'],fitTable['b'],fitTable['d'],fitTable['c']).copy()
    IC50valLum=IC50valLum.reset_index()
    #This is the variable containing the IC50 Value for the cell line and drug combo
    IC50valLum2=IC50valLum.iat[0,1]
    IC50df = IC50df.append({'Cell_Line':CL,'Drug':dd,'IC50_Lum':IC50valLum2}, ignore_index=True)
    IC50lum =r'IC50 : {:0.3e}'.format(IC50valLum2)+' M'
    refDose = -np.linspace(min(full_df['DClog'])*1.2,max(full_df['DClog'])*0.8,256)
    refDose = (10**-refDose)
    sns.lmplot(x='DClog',y='DIMrNorm',data=full_df,fit_reg=False)        
    for fit in fitData:
        plt.plot([pDose(i) for i in refDose],[ll4(i,*[fit[i] for i in ['b','c','d','e']]) for i in refDose])
    sns.set(style='darkgrid')            
    plt.ylabel('Normalized Luminescence Rate')
    plt.xlabel('log Drug Concentration (M)')
    plt.title(CL+' treated with '+dd)
    plt.legend(title=None, loc='upper right', labels=[IC50lum])
    plt.savefig('/Users/claytonwandishin/December 14 RT glow run/plots/PaperFigs/DRC_Lum/H1048/AlgoV2/DRC/'+CL+'_'+dd+'.png', dpi=300, bbox_inches='tight')
    plt.show()
full_df = full_df.merge(IC50df)
#####THIS GENERATES THE R2 vs. Slice size dataframe
r2df = pd.DataFrame()
r2df['DrugCon'] = dcrrlabel
r2df['Well'] = dcrrwell
r2df['SliceSize'] = dcrrslice
r2df['R2_Value'] = dcrr

well_list_r2df = r2df['Well'].unique()

def set_custom_palette(series, max_color = 'tab:orange', other_color = 'tab:blue'):
    max_val = pd.Series(series).max()
    print(max_val)
    pal = []
    
    for item in series:
        if item == max_val:
            pal.append(max_color)
        else:
            pal.append(other_color)
        print(pal)
    return pal

def set_custom_markers(series, max_mark = '^', other_mark = 'o'):
    max_val = series.max()
    mar = []
    
    for item in series:
        if item == max_val:
            mar.append(max_mark)
        else:
            mar.append(other_mark)
    return mar

def set_custom_marker_size(series, max_mar_size = 200, other_mar_size = 50):
    max_val = series.max()
    marsize = []
    
    for item in series:
        if item == max_val:
            marsize.append(max_mar_size)
        else:
            marsize.append(other_mar_size)
    return marsize

for uw in well_list_r2df:
    r2dfplot = r2df.loc[r2df['Well'] == uw]
    cuspal=set_custom_palette(r2dfplot['R2_Value'])
    cusmar=set_custom_markers(r2dfplot['R2_Value'])
    cussize=set_custom_marker_size(r2dfplot['R2_Value'])
    r2dfplot['cpallette']=cuspal
    sns.set(style='darkgrid')
    sns.scatterplot(data=r2dfplot, x=r2dfplot['SliceSize'], y=r2dfplot['R2_Value'], hue=r2dfplot['R2_Value'],hue_order=r2dfplot['R2_Value'], palette=cuspal,style=r2dfplot['R2_Value'],markers=cusmar, style_order=r2dfplot['R2_Value'], size = r2dfplot['R2_Value'], sizes=cussize,size_order=r2dfplot['R2_Value'], legend=False)
    plt.ylabel('')
    plt.xlabel('')
    plt.xlim(19,3)
    #plt.title(CL+' well '+uw)
    #plt.savefig('/Users/claytonwandishin/December 14 RT glow run/plots/PaperFigs/DRC_Lum/H1048/AlgoV2/'+CL+'_'+dd+'/'+CL+'_well_'+uw+'R2plot.png', dpi=300, bbox_inches='tight')
    plt.show()
#pmdf = pd.read_csv('/Users/claytonwandishin/December 14 RT glow run/20201216_Lum_CellCounts_TOTAL_DIPcomp.csv')
#full_df = full_df.merge(pmdf)


#full_df.to_csv('/Users/claytonwandishin/December 14 RT glow run/20201216_Lum_CellCounts_TOTAL_DIM_33.csv')
'''
dclist_r2df = r2df['DrugCon'].unique()
for d in dclist_r2df:
    r2dfplot = r2df.loc[r2df['DrugCon'] == d]
    sns.set(style='darkgrid')
    sns.scatterplot(data=r2dfplot, x=r2dfplot['SliceSize'], y=r2dfplot['R2_Value'])
    plt.xlim(19,3)
    plt.ylabel('$R^2 Value')
    plt.xlabel('Timepoint Slice Size (points)')
    plt.title(CL+' with '+str(d)+' Barasertib')
    plt.show()
'''

######################################################
######################################################
#######################################################
#######################################################

# This part creates rate-based DataFrames so that the DRC generated from Lum and Direct Counting can be compared
'''
NoDrugLumandCounts = full_df.loc[full_df['Drug_Conc'] == 0]
NDLCC = NoDrugLumandCounts
NDLCC['CC+1']=NDLCC['Live_Dead']+1
NDLCC['log2_Live_Dead+1']=np.log2(NDLCC['CC+1'])

Lnorm = NDLCC.loc[:,'log2_Live_Dead+1']
Lnorm = pd.DataFrame(Lnorm)
Lnorm['TotHour_Image']= NDLCC.loc[:,'TotHour_Image']

DrugCon = ['0','0.00000996','0.00000249','0.000000623','0.000000156','0.0000000389','0.00000000934','0.00000000233','0.000000000623','0.000000000156','0.0000000000778']

DrugConfloat = [float(x) for x in DrugCon]

cell_line_list = ['H1048']
#drug_list_fig = ['AZD-1152','SNS-314','trametinib']
drug_list_fig = ['barasertib']
#drug_list_fig = completeDF['Drug'].unique()
#etoposide and DMSO fail
#remove 0 in Lum for SCH
#remove 0 for vemurafinib
#remove 7.78e-11 Lum and 9.96 e-06 for direct
#remove 2.49e-6 and 7.78e-11 for azd1152
DrugConfloat.sort()
DrugConfloatDIRECT = DrugConfloat.copy()
#DrugConfloatDIRECT.remove(0)
#DrugConfloatDIRECT.remove(2.33e-09)
#DrugConfloat.remove(2.49e-06)
#DrugConfloat.remove(0)
#DrugConfloatDIRECT.remove(1.56e-7)
#DrugConfloat.remove(9.96e-06)
#DrugConfloat.remove(7.78e-11)
DrugConfloatDIRECT.sort()
DipDF = pd.DataFrame()
DipDF['Drug_Con']=[]
DipDF['DIP_Rate']=[]
DipDF['Count_Type']=[]
DipDF['DIP_Rate_Norm']=[]
DipDF['Response_Ratio']=[]
DipDFw = []
count_type = ['Direct','Lum']
for CLL in cell_line_list:
    for DD in drug_list_fig:
        oneCLdf = full_df.loc[full_df['Cell_Line'] == CLL]
        OCLDF = oneCLdf
        onedrugdf = OCLDF.loc[OCLDF['Drug'] == DD]
        ODDF = onedrugdf
        ODDFt100 = ODDF.loc[ODDF['TotHour_Image']<maxtime0]
        ODDftLum = ODDF.loc[ODDF['TotHour_Lum']<maxtime0]
       

        #DrugConfloat = completeDF['Drug_Conc'].unique()
        
        for c in DrugConfloatDIRECT:
            indivdrugcondip = ODDFt100.loc[ODDFt100['Drug_Conc'] == c]
            slope, intercept, r_value, pv, se =     stats.linregress(indivdrugcondip['TotHour_Image'],indivdrugcondip['log2_Live_Dead'])
            DipDF=DipDF.append([{'Drug_Con':c, 'DIP_Rate':slope, 'Count_Type':'Direct'}], ignore_index=True)
            print(str(c)+" "+str(slope)+"Direct") 
            
            
            
        for c in DrugConfloat:
            indivdrugcondip = ODDftLum.loc[ODDftLum['Drug_Conc'] == c]
            slope, intercept, r_value, pv, se =     stats.linregress(indivdrugcondip['TotHour_Lum'],indivdrugcondip['RLU'])
            DipDF=DipDF.append([{'Drug_Con':c, 'DIP_Rate':slope, 'Count_Type':'Lum'}], ignore_index=True)
            print(str(c)+" "+str(slope)+"Lum") 




#This needs to be fixed to reflect the different scales of lum and direct counting
            LumDip = DipDF.loc[DipDF['Count_Type'] == 'Lum']
        
            
            DirectDip = DipDF.loc[DipDF['Count_Type'] == 'Direct']
        
            LumDipMax = LumDip['DIP_Rate'].max()
            DirectDipMax = DirectDip['DIP_Rate'].max()
            for i, row in DipDF.iterrows():
                if DipDF.at[i,'Count_Type'] == 'Lum':
                    DipNorm = DipDF.at[i,'DIP_Rate']/LumDipMax
                    DipDF.at[i,'DIP_Rate_Norm'] = DipNorm
                    
                    
            for i, row in DipDF.iterrows():
                if DipDF.at[i,'Count_Type'] == 'Direct':
                    DipNorm = DipDF.at[i, 'DIP_Rate']/DirectDipMax
                    DipDF.at[i,'DIP_Rate_Norm'] = DipNorm
            
            LumDipNormMin = LumDip['DIP_Rate_Norm'].min()
            LumDipNormMax = LumDip['DIP_Rate_Norm'].max() 
            DirectDipNormMin = DirectDip['DIP_Rate_Norm'].min()
            DirectDipNormMax = DirectDip['DIP_Rate_Norm'].max()  
            for i, row in DipDF.iterrows():
                if DipDF.at[i,'Count_Type'] == 'Lum':
                    DipRR = (DipDF.at[i,'DIP_Rate_Norm']- LumDipNormMin)/(LumDipNormMax - LumDipNormMin)
                    DipDF.at[i,'Response_Ratio'] = DipRR
                    
                    
            for i, row in DipDF.iterrows():
                if DipDF.at[i,'Count_Type'] == 'Direct':
                    DipRR = (DipDF.at[i,'DIP_Rate_Norm']- DirectDipNormMin)/(DirectDipNormMax - DirectDipNormMin)
                    DipDF.at[i,'Response_Ratio'] = DipRR
            
            
        DipDF['DClog']=pDose(DipDF['Drug_Con']+.0000000001)
        
        #DipDF['Response_Ratio']= DipDF.groupby(['Count_Type']).transform(lambda x: (DIP_Rate_Norm - DIP_Rate_Norm.min())/ LumDipDIP_Rate_Norm.max() - DIP_Rate_Norm.min()))
#DipDF['DClog']=np.log(DipDF['Drug_Con']+.0000000001)


#FITTING FUNCTION NEEDS TO BE SPLIT INTO TWO
            #DipDFw=DipDF.loc[DipDF['Count_Type'] == f]
        compoundData = DipDF.groupby(['Count_Type'])
        fitData = []
        for name,group in compoundData:
            fitCoefs, covMatrix = opt.curve_fit(ll4, group.Drug_Con,         group.Response_Ratio)
            resids = group.Response_Ratio-group.Drug_Con.apply(lambda x: ll4(x,*fitCoefs))
            curFit = dict(zip(['b','c','d','e'],fitCoefs))
            curFit['Count_Type']=name
            curFit['residuals']=sum(resids**2)
            fitData.append(curFit)
            fitCompound = [ item['Count_Type'] for item in fitData]
            fitTable = pd.DataFrame(fitData).set_index('Count_Type')
            
            if name == 'Direct':
                IC50valDirect =       IC50(fitTable['e'],fitTable['b'],fitTable['d'],fitTable['c']).copy()
                IC50valDirect=IC50valDirect.reset_index()
                IC50valDirect2=IC50valDirect.at[0,0]
                IC50direct = r'IC50 : {:0.3e}'.format(IC50valDirect2)+' M'
            if name == 'Lum':
                IC50valLum =       IC50(fitTable['e'],fitTable['b'],fitTable['d'],fitTable['c']).copy()
                IC50valLum=IC50valLum.reset_index()
                IC50valLum2=IC50valLum.at[1,0]
                IC50lum =r'IC50 : {:0.3e}'.format(IC50valLum2)+' M'
            refDose = -np.linspace(min(DipDF['DClog'])*0.6,max(DipDF['DClog'])*1.5,256)
            refDose = (10**-refDose)*1e6
            sns.lmplot(x='DClog',y='Response_Ratio',data=DipDF,hue='Count_Type',fit_reg=False)        
            for fit in fitData:
                plt.plot([pDose(i) for i in refDose],[ll4(i,*[fit[i] for i in ['b','c','d','e']]) for i in refDose])



                
#IC50lum =r'IC50 : {:0.3e}'.format(IC50valLum2.iloc[0])+' M'
#IC50direct = r'IC50 : {:0.3e}'.format(IC50valDirect2.iloc[0])+' M'
sns.set(style='darkgrid')            
plt.ylabel('Response Ratio')
plt.xlabel('log Drug Concentration (M)')
plt.title(CL+" treated with "+DD)
plt.legend(title=None, loc='upper right', labels=[IC50direct,IC50lum])
#plt.show()
#plt.savefig('/Users/claytonwandishin/December 14 RT glow run/plots/PaperFigs/DRC_Lum/'+CL+'/'+CL+'DirectandLumCOMPARISON_DRC_RESPONSERATIO_RAWlum'+DD+'DECrun.png', dpi=300, bbox_inches='tight')
#ODDF.to_csv('/Users/claytonwandishin/December 14 RT glow run/plots/PaperFigs/DRC_Lum/'+CL+'/'+CL+'YM155_OneDrugDF.csv')
'''