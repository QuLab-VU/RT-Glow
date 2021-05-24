from pysb import *
from pysb.simulator import ScipyOdeSimulator
import numpy as np
import matplotlib.pyplot as plt
from networkx.algorithms.bipartite.tests.test_matrix import sp

Model()

# Monomer('Cell')
# Parameter('Cell_0', 100)
# Initial(Cell(), Cell_0)
# Parameter('k_div', 1)
# Rule('Cell_division', Cell() >> Cell() + Cell(), k_div)
# Parameter('k_dth', 0.9)
# Rule('Cell_death', Cell() >> None, k_dth)
# Observable('Cell_tot', Cell())
#  
# tspan = np.linspace(0,100,1001)
# sim = ScipyOdeSimulator(model, tspan, verbose=True)
# x = sim.run()
#  
# plt.plot(tspan, np.log(x.observables['Cell_tot']), 'k--', lw=2, label='Total')
# plt.legend(loc=0)
# plt.show()
#  
# quit()

Monomer('Cell', ['c', 'c'])

Parameter('Cell_0', 100)
Initial(Cell(c=MultiState(None,None)), Cell_0)

Parameter('kf_bind', 1)
Parameter('kr_bind', 1)
Rule('Cell_Cell_bind', Cell(c=None) + Cell(c=MultiState(None,None)) | 
     Cell(c=1) % Cell(c=MultiState(1,None)), kf_bind, kr_bind)

Parameter('k_div', 1)
Rule('Cell_division', Cell(c=None) >> Cell(c=1) % Cell(c=MultiState(1,None)), k_div)

Parameter('k_dth', 0.9)
Rule('Cell_death', Cell(c=MultiState(None,None)) >> None, k_dth)

Observable('Cell_tot', Cell())

from pysb.bng import generate_equations
generate_equations(model, verbose=True, max_agg=10)

tspan = np.linspace(0,100,1001)
sim = ScipyOdeSimulator(model, tspan, verbose=True)
x = sim.run()

plt.figure()
for i,sp in enumerate(model.species):
    plt.plot(tspan, x.all['__s%d' % i], label='size %d' % (i+1), lw=2)
plt.plot(tspan, x.observables['Cell_tot'], 'k--', lw=2, label='Total')
exp_growth = np.array(Cell_0.value*np.exp((k_div.value-k_dth.value)*tspan))
plt.plot(tspan, exp_growth, 'r--', lw=2, label=r'y=%de$^{%.1ft}$' % 
         (Cell_0.value, k_div.value-k_dth.value))
plt.xlabel('Time')
plt.ylabel('Cell count')
plt.ylim(ymin=0, ymax=25000)
plt.legend(loc=0)

plt.figure()
for i,sp in enumerate(model.species):
    plt.plot(tspan, np.log(x.all['__s%d' % i]), label='size %d' % (i+1), lw=2)
plt.plot(tspan, np.log(x.observables['Cell_tot']), 'k--', lw=2, label='Total')
plt.plot(tspan, np.log(exp_growth), 'r--', lw=2, label=r'ln y=ln %d+%.1ft' % 
         (Cell_0.value, k_div.value-k_dth.value))
plt.xlabel('Time')
plt.ylabel('ln(Cell count)')
plt.legend(loc=0, ncol=3)

plt.show()

# for sp in model.species:
#     print(sp)
# print()
# for rxn in model.reactions:
#     print(rxn)

# print('nspecies', len(model.species))
# print('nrxns', len(model.reactions))
