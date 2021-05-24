from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import numpy as np

plt.figure('DIP rate dose-response curve')
plt.title('NCIH1048 + TAK-901 DIP dose-response curve')

# SIMULATED DATA

sim_data = np.genfromtxt('NCIH1048_TAK-901/simulated_DIP_rates_Thunor.csv', 
                     delimiter=',', names=True) #, encoding)
print(sim_data)
print(sim_data.dtype.names)

sim_control = [d['DIP'] for d in sim_data if d['Dose'] < 2e-11]
sim_control_mean = np.mean(sim_control)

plt.plot(np.log10(sim_data['Dose']), sim_data['DIP'], #/sim_control_mean, 
         "r^", mfc='None', label='sim data')

def Hill_direct_effect(drug, Emax, E0, EC50, h):
#     E0 = 1
    return Emax + (E0-Emax)/(1+(drug/EC50)**h)

popt, pcov = curve_fit(Hill_direct_effect, sim_data['Dose'], sim_data['DIP'], #/sim_control_mean, 
                       maxfev = 10000000)

conc_points = 10**np.linspace(-5, -11, 61)
sim_fit = [Hill_direct_effect(d, popt[0], popt[1], popt[2], popt[3]) for d in conc_points]
print('sim fit: Emax = %.2g, E0 = %.2g, EC50 = %.2g, h = %.2g' 
      % (popt[0], popt[1], popt[2], popt[3]))

# IC50 = (0.5/(0.5-Emax))^(1/h)*EC50
# sim_IC50 = (0.5/(0.5-popt[0]))**(1./popt[3])*popt[2]
sim_IC50 = (0.5*popt[1]/(0.5*popt[1]-popt[0]))**(1./popt[3])*popt[2]

plt.plot(np.log10(conc_points), sim_fit, 'r', label='sim IC50 = %.2g' % sim_IC50)
plt.xlabel('log10(dose M)')
plt.ylabel('DIP rate (1/h)')

# EXPERIMENTAL DATA

exp_data = np.genfromtxt('experimental_DIP_rates_Thunor.csv', 
                     delimiter=',', names=True) 
print(exp_data)
print(exp_data.dtype.names)

exp_control = [d['DIP'] for d in exp_data if d['Dose'] < 2e-11]
exp_control_mean = np.mean(exp_control)

plt.plot(np.log10(exp_data['Dose']), exp_data['DIP'], #/exp_control_mean, 
         "ks", mfc='None', label='exp data')

popt, pcov = curve_fit(Hill_direct_effect, exp_data['Dose'], exp_data['DIP'], #/exp_control_mean, 
                       maxfev = 10000000)

exp_fit = [Hill_direct_effect(d, popt[0], popt[1], popt[2], popt[3]) for d in conc_points]

# exp_IC50 = (0.5/(0.5-popt[0]))**(1./popt[3])*popt[2]
exp_IC50 = (0.5*popt[1]/(0.5*popt[1]-popt[0]))**(1./popt[3])*popt[2]

plt.plot(np.log10(conc_points), exp_fit, 'k', label='exp IC50 = %.2g' % exp_IC50)
plt.xlabel('log10(dose M)')
plt.ylabel('DIP rate (1/h)')

plt.legend(loc=0)
plt.show()

