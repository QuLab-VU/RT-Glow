from __future__ import print_function
from pysb import *
import numpy as np
from scipy.constants import N_A

Model()

Monomer('Cells', ['NADPHoxidoreductase1', 'drug'], {'drug': ['drugged', 'undrugged']})
Monomer('Drug')

Monomer('NanoLuc', ['BindingSite'])
Monomer('ProSubstrate', ['BindingSiteNADPH_OR_and_NL', 'state'], {'state': ['ProSub', 'Substrate']})


#'''
#starting parameters
Parameter('K3', 2.011170567427241320e-04)   # estimated
Parameter('K5', 1.0) # [1/day] ... time unit!  leave
#Parameter('K6', 1-np.log(2)/(39.2/24)) #DIP-rate = 1/39.2h number taken out of the literature
Parameter('K6', 1-np.log(2)/(48.46/24)) #doubling time taken out of the actual cell count data   leave
Parameter('Km', 5.798362378769308270e+04)   # estimated
Parameter('Kcat', 6.398140950482379594e+03)   # estimated
Parameter('cc', 6500) # define the carrying capacity to include confluency    leave

Parameter('kon', 1000.0)     # estimated
Parameter('koff', 1000.0)    # estimated
Parameter('K5d', 1.0)     # estimated -> encode in costfunction (if check fails - has to be thrown out) k5d <= k5
Parameter('K6d', 1.0)     # estimated k6d >= k6
#additional constraint: k5d-k6d <= k5-k6 (proliferation rate)

#THIS PARAMETER WILL NOT BE FIT (THE OTHERS WILL)
Parameter('Gain', 8.0)


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
#Expression('confluency', (K5-K6)/cc)
Observable('Live_Cells', Cells())
Expression('confluency', Live_Cells*(K5-K6)/cc)
#Rule('carrying_capacity', Cells()+Cells() >> Cells(), confluency)
Rule('carrying_capacity', Cells() >> None, confluency)



Parameter('Cells_0', 1572)
Parameter('NanoLuc_0', 209424)
#Parameter('NanoLuc_0', 1.26073248E16)
Parameter('ProSubstrate_0', 9e8)
Parameter('Drug_0', 1.0e-6*N_A*80e-6)   #1.0e-6 =1uM -> different drug concentration via multiple objective

Initial(Cells(NADPHoxidoreductase1=None, drug='undrugged'), Cells_0)
Initial(Drug(), Drug_0)
print(Drug(), type(Drug()))
Initial(NanoLuc(BindingSite=None), NanoLuc_0)
Initial(ProSubstrate(BindingSiteNADPH_OR_and_NL=None, state='ProSub'), ProSubstrate_0)



Observable('Unconverted_Prosubstrate', ProSubstrate(BindingSiteNADPH_OR_and_NL=None, state='ProSub'))

Observable('Converted_Prosubstrate', ProSubstrate(BindingSiteNADPH_OR_and_NL=None, state='Substrate'))



'''
Monomer('A')
Parameter('A0', 1000)
Parameter('k', 1)

Initial(A(), A0)

Rule('Adeg', A() >> None, k)


x = sim.run(tspan = tspan)  # -> A = 10

y = sim.run(tspan = tspan, param_values = [100, 2], initials = x.species[-1])'''