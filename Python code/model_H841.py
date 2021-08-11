from __future__ import print_function
from pysb import *
import numpy as np
from scipy.constants import N_A

Model()

Monomer('Cells', ['NADPHoxidoreductase1', 'drug'], {'drug': ['drugged', 'undrugged']})
Monomer('Drug')

Monomer('NanoLuc', ['BindingSite'])
Monomer('ProSubstrate', ['BindingSiteNADPH_OR_and_NL', 'state'], {'state': ['ProSub', 'Substrate']})

#starting parameters
Parameter('K3', 2.011170567427241320e-04)   # estimated
Parameter('K5', 1.0) # [1/day] ... time unit!  not estimated
Parameter('K6', 1-np.log(2)/(48.46/24)) # doubling time taken out of the actual cell count data -> not estimated
Parameter('Km', 5.798362378769308270e+04)   # estimated
Parameter('Kcat', 6.398140950482379594e+03)   # estimated
Parameter('cc', 6500) # define the carrying capacity to include confluency -> not estimated

Parameter('kon', 1000.0)     # estimated
Parameter('koff', 1000.0)    # estimated
Parameter('K5d', 1.0)     # estimated -> encode in costfunction (if check fails - has to be thrown out) k5d <= k5
Parameter('K6d', 1.0)     # estimated k6d >= k6

Parameter('Gain', 8.0)   # not estimated

#=====drug reaction for both cell states (w/wo NADPHoxidoreductase1) react with drug
Rule('DrugReaction', Cells(drug='undrugged') + Drug() >> Cells(drug='drugged') + Drug(), kon)
Rule('DrugOffCell', Cells(drug='drugged') >> Cells(drug='undrugged'), koff)
Rule('CellDivisionNoDrug', Cells(drug='undrugged') >> Cells(drug='undrugged') + Cells(NADPHoxidoreductase1=None, drug='undrugged'), K5)
Rule('CellDivisionDrug', Cells(drug='drugged') >> Cells(drug='drugged') + Cells(NADPHoxidoreductase1=None, drug='drugged'), K5d)
Rule('CellDeathNoDrug', Cells(NADPHoxidoreductase1=None, drug='undrugged') >> None, K6)
Rule('CellDeathDrug', Cells(NADPHoxidoreductase1=None, drug='drugged') >> None, K6d)

Rule('ProSubtoSubConversion',Cells(NADPHoxidoreductase1=None) + ProSubstrate(BindingSiteNADPH_OR_and_NL=None, state='ProSub', ) >> Cells(NADPHoxidoreductase1=None) + ProSubstrate(BindingSiteNADPH_OR_and_NL=None, state='Substrate'), K3)
Observable('Substrate', ProSubstrate(BindingSiteNADPH_OR_and_NL=None, state='Substrate'))
Expression('NL_MM', Kcat / (Km + Substrate))
Rule('NanolucBinding', NanoLuc(BindingSite=None) + ProSubstrate(BindingSiteNADPH_OR_and_NL=None, state='Substrate') >> NanoLuc(BindingSite=None), NL_MM)
Expression('Luminescence', Gain * Substrate)
Observable('Live_Cells', Cells())
Observable('Undrugged_Live_Cells', Cells(drug='undrugged'))
Observable('Drugged_Live_Cells', Cells(drug='drugged'))

Expression('confluency_nodrug', Live_Cells*(K5-K6)/cc)
Rule('carrying_capacity_nodrug', Cells(drug='undrugged') >> None, confluency_nodrug)
Expression('confluency_drug', Live_Cells*(K5d-K6d)/cc)
Rule('carrying_capacity_drug', Cells(drug='drugged') >> None, confluency_drug)

Parameter('Cells_0', 1572)
Parameter('NanoLuc_0', 209424)
Parameter('ProSubstrate_0', 9e8)
Parameter('Drug_0', 0)

Initial(Cells(NADPHoxidoreductase1=None, drug='undrugged'), Cells_0)
Initial(Drug(), Drug_0)
Initial(NanoLuc(BindingSite=None), NanoLuc_0)
Initial(ProSubstrate(BindingSiteNADPH_OR_and_NL=None, state='ProSub'), ProSubstrate_0)

Observable('Unconverted_Prosubstrate', ProSubstrate(BindingSiteNADPH_OR_and_NL=None, state='ProSub'))
Observable('Converted_Prosubstrate', ProSubstrate(BindingSiteNADPH_OR_and_NL=None, state='Substrate'))

'''
############
from pysb.simulator import ScipyOdeSimulator
import matplotlib.pyplot as plt

data = np.genfromtxt('H841_TAK901_Lum.csv', delimiter=',', names=True)
print(data.dtype.names)

rate_mask = np.array([True, False, False, True, True, False, True, True, True, True, False, False, False, False, False])
rate_mask = [i for i,r in enumerate(rate_mask) if r]

pvals = np.genfromtxt('gathered_params_H841.txt', delimiter=',')

bkgrd = 118.6666667
for name in data.dtype.names[1:]:
    data[name] -= bkgrd
    data[name] = np.clip(data[name], a_min=0, a_max=None)

sim = ScipyOdeSimulator(model, verbose=True)

param_values = {}
for i in range(len(rate_mask)):
    param_values[model.parameters[rate_mask[i]].name] = pvals[-1][i]

# for key in param_values.keys():
#     print(key, param_values[key])

# K3 7.199958902647853e-05
# Km 12549465.660542943
# Kcat 5097153.474404186

# Higher K3 -> Higher lum
# Higher Km -> Higher lum
# Higher Kcat -> Lower lum

### MANUALLY ADJUST PARAMETERS
param_values['K3'] = 7e-6
param_values['Km'] = 5e9
param_values['Kcat'] = 5e7
##############################

## PRE-SIMULATION
tspan = [0, 0.5]
x = sim.run(tspan=tspan, 
            param_values=param_values,
            initials = {Cells(NADPHoxidoreductase1=None, drug='undrugged') : 300,
                        Drug() : 0})

# REGULAR SIMULATION
tspan = np.linspace(0.5, 7, 501)
initials = x.species[-1]
x = sim.run(tspan=tspan, 
            param_values=param_values,
            initials = initials)

plt.figure()
plt.plot(tspan, x.observables['Live_Cells'], label='Live Cells')
plt.xlabel('Time (day)')
plt.ylabel('Count')
plt.legend(loc=0)

plt.figure()
for name in data.dtype.names[1:]:
    plt.plot(data['TotHour'][1:]/24, data[name][1:], '.', label=name)
plt.plot(tspan, x.expressions['Luminescence'], label='Simulation')
plt.xlabel('Time (day)')
plt.ylabel('Luminescence')
plt.legend(loc=0)

plt.show()
'''
