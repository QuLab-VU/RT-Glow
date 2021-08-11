from pysb import *
from pysb.simulator import ScipyOdeSimulator
import numpy as np
import matplotlib.pyplot as plt

Model()

# Cell + Drug -> Cell* + Drug  kf
# Cell* -> Cell  kr

# Cell -> 2 Cell    kdiv_u
# Cell -> None      kdth_u
# 2 Cell -> Cell    (kdiv_u-kdth_u)/CC

# Cell* -> 2 Cell*  kdiv_d
# Cell* -> None     kdth_d
# 2 Cell* -> Cell*  (kdiv_u-kdth_u)/CC

# Pro + Cell -> Sub + Cell  k_pro_to_sub
# Sub + E -> E  kcat/(Km + [Sub])

Monomer('Cell', ['state'], {'state' : ['u', 'd']})
Monomer('Drug')
Monomer('Pro')
Monomer('Sub')
Monomer('Enzyme')

Initial(Cell(state='u'), Parameter('Cell_0', 100))
Initial(Drug(), Parameter('Drug_0', 0))
Initial(Pro(), Parameter('Pro_0', 1e6)) # ???

Observable('Cell_tot', Cell())
Observable('Light_tot', Light())
Observable('Pro_tot', Pro())
Observable('Sub_tot', Sub())

Parameter('kf', 1)
Parameter('kr', 1)

Rule('Cell_drug_forward', Cell(state='u') + Drug() >> Cell(state='d') + Drug(), kf)
Rule('Cell_drug_reverse', Cell(state='d') >> Cell(state='u'), kr)

Parameter('cc', 10000) # carrying capacity

Parameter('k_div_u', 1) # /day 0.3-1
Parameter('k_dth_u', 0)
Rule('Cell_div_u', Cell(state='u') >> Cell(state='u') + Cell(state='u'), k_div_u)
Rule('Cell_dth_u', Cell(state='u') >> None, k_dth_u)

Parameter('k_div_d', 1) 
Parameter('k_dth_d', 1.1)
Rule('Cell_div_d', Cell(state='d') >> Cell(state='d') + Cell(state='d'), k_div_d)
Rule('Cell_dth_d', Cell(state='d') >> None, k_dth_d)

Parameter('k_crowd', (k_div_u.value-k_dth_u.value)/cc.value)
Rule('Cell_crowd', Cell() + Cell() >> Cell(), k_crowd)

Parameter('k_pro_to_sub', 1e-5)
Rule('Pro_to_Sub', Pro() + Cell() >> Sub() + Cell(), k_pro_to_sub)

Parameter('kcat_Et', 1e5)
Parameter('Km', 1)
Expression('rate_emit', kcat_Et/(Km + Sub_tot))
Rule('Sub_emits_Light', Sub() >> None, rate_emit)

tspan = np.linspace(0, 3, 301)
sim = ScipyOdeSimulator(model, verbose=True)

#####
# x = sim.run(tspan, param_values={'k_div_u' : 0, 'k_crowd_u' : 0}).observables
#  
# plt.plot(tspan, x['Pro_tot'], lw=2, label='Pro')
# plt.plot(tspan, x['Sub_tot'], lw=2, label='Sub')
# plt.plot(tspan, x['Light_tot'], lw=2, label='Light')
# plt.xlabel('time (d)')
# plt.ylabel('amount')
# plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
#  
# plt.legend(loc=0)
# plt.tight_layout()
# plt.show()
#  
# quit()

#####
Cell_tot = np.linspace(1e2,1e4,10,endpoint=True)
Light_tot = []
plt.figure()
for cell_0 in Cell_tot:
    print(cell_0)
    x = sim.run(tspan, initials={Cell(state='u') : cell_0}, param_values={'k_div_u' : 0, 'k_crowd' : 0}).observables
    plt.plot(tspan, x['Light_tot'], lw=2, label="%g" % cell_0)
    Light_tot.append(x['Light_tot'][-1])
plt.xlabel('time (d)')
plt.ylabel('light')
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.legend(loc=0)
plt.tight_layout()
     
plt.figure()
plt.plot(Cell_tot, Light_tot, 'ro')
plt.xlabel('total cells')
plt.ylabel('light')
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.tight_layout()
 
# plt.show()
# quit()
#####

x = sim.run(tspan=np.linspace(0, 10, 101))

fig = plt.figure()
ax1 = fig.add_subplot(111)
# line1 = ax1.plot(tspan, np.log2(x['Cell_tot']/x['Cell_tot'][0]), 'b', lw=2, label='Cells')
# ax1.set_ylabel('doublings')
lines = ax1.plot(x.tout[0], x.observables['Cell_tot'], 'b', lw=2, label='Cells')
ax1.set_ylabel('count')
ax1.set_xlabel('time (day)')
ax1.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

ax2 = fig.add_subplot(111, sharex=ax1, frameon=False)
lines += ax2.plot(x.tout[0], x.observables['Light_tot'], 'r', lw=2, label='Light')
ax2.set_ylabel('light')
ax2.yaxis.tick_right()
ax2.yaxis.set_label_position("right")
ax2.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

#### In silico "experimental" data
from scipy.stats import norm
t_expt = [x.tout[0][i] for i in range(0,len(x.tout[0]),5)]
light_expt = [x.observables['Light_tot'][i] for i in range(0,len(x.tout[0]),5)]
light_expt = [norm.rvs(l,0.075*l) for l in light_expt]
lines += ax2.plot(t_expt, light_expt, 'kx', mfc = "none", mew=2, label='"experiment"')
####

# ax1.legend(loc=0)
labels = [l.get_label() for l in lines]
ax1.legend(lines, labels, loc=0)

fig.tight_layout()

plt.show()
