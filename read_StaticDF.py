import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import Normalize
from scipy.stats import linregress
from scipy.optimize import curve_fit

infile = 'data/20201216_Lum_Image_Run/StaticDF.csv'

raw_data = np.genfromtxt(infile, names=True, dtype=None, delimiter=',', encoding='UTF-8')
# print(raw_data)
print(raw_data.dtype.names)

cell_lines = np.unique([d['Cell_Line'] for d in raw_data])

print(cell_lines)

#####
cell_line = 'H1048'
#####

exp_data = np.array([d for d in raw_data if d['Cell_Line'] == cell_line])

def linear_zero_intercept(x, m):
    return m*x

def linear_unit_slope(x, b):
    return x + b

plt.figure()
plt.plot(exp_data['Cell_Conc'], exp_data['RLU'], 'o', mfc='None')
fit1 = linregress(exp_data['Cell_Conc'], exp_data['RLU'])
label = 'y = %g*x + %g (R$^2$ = %g)' % (fit1.slope, fit1.intercept, fit1.rvalue**2)
plt.plot(exp_data['Cell_Conc'], fit1.slope * exp_data['Cell_Conc'] + fit1.intercept,
         lw=2, label=label)
popt, pcov = curve_fit(linear_zero_intercept, exp_data['Cell_Conc'], exp_data['RLU'])
# calculate R^2 #####
residuals = exp_data['RLU'] - linear_zero_intercept(exp_data['Cell_Conc'], *popt)
ss_res = np.sum(residuals**2)
ss_tot = np.sum((exp_data['RLU']-np.mean(exp_data['RLU']))**2)
r_squared = 1 - (ss_res / ss_tot)
#####################
label = 'y = %g*x (R$^2$ = %g)' % (popt[0], r_squared)
plt.plot(exp_data['Cell_Conc'], popt[0] * exp_data['Cell_Conc'], lw=2, label=label)
plt.xlabel('cell count')
plt.ylabel('luminescence (RLU)')
plt.legend(loc=0)
plt.tight_layout()

plt.figure()
plt.plot(np.log2(exp_data['Cell_Conc']), np.log2(exp_data['RLU']), 'o', mfc='None')
fit2 = linregress(np.log2(exp_data['Cell_Conc']), np.log2(exp_data['RLU']))
label = 'y = %g*x + %g (R$^2$ = %g)' % (fit2.slope, fit2.intercept, fit2.rvalue**2)
plt.plot(np.log2(exp_data['Cell_Conc']), fit2.slope * np.log2(exp_data['Cell_Conc']) + fit2.intercept,
         lw=2, label=label)
popt, pcov = curve_fit(linear_unit_slope, np.log2(exp_data['Cell_Conc']), np.log2(exp_data['RLU']))
# calculate R^2 #####
residuals = np.log2(exp_data['RLU']) - linear_unit_slope(np.log2(exp_data['Cell_Conc']), *popt)
ss_res = np.sum(residuals**2)
ss_tot = np.sum((np.log2(exp_data['RLU'])-np.mean(np.log2(exp_data['RLU'])))**2)
r_squared = 1 - (ss_res / ss_tot)
#####################
label = r'y = x + %g (R$^2$ = %g)' % (popt[0], r_squared)
plt.plot(np.log2(exp_data['Cell_Conc']), popt[0] + np.log2(exp_data['Cell_Conc']), lw=2, label=label)
plt.xlabel('log2(cell count)')
plt.ylabel('log2(luminescence (RLU))')
plt.legend(loc=0)
plt.tight_layout()

plt.show()
