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
Parameter('K3d', 2.011170567427241320e-04)   # estimated
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

Rule('ProSubtoSubConversionNoDrug',Cells(drug='undrugged', NADPHoxidoreductase1=None) + ProSubstrate(BindingSiteNADPH_OR_and_NL=None, state='ProSub', ) >> Cells(drug='undrugged', NADPHoxidoreductase1=None) + ProSubstrate(BindingSiteNADPH_OR_and_NL=None, state='Substrate'), K3)
Rule('ProSubtoSubConversionDrug',Cells(drug='drugged', NADPHoxidoreductase1=None) + ProSubstrate(BindingSiteNADPH_OR_and_NL=None, state='ProSub', ) >> Cells(drug='drugged', NADPHoxidoreductase1=None) + ProSubstrate(BindingSiteNADPH_OR_and_NL=None, state='Substrate'), K3d)
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
