from model_H841 import model as m841
from pysb.simulator import ScipyOdeSimulator
import csv
import math as m
from scipy import stats
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import N_A
from statistics import mean


cell_ct1 = []
cell_ct2 = []
cell_ct3 = []
cell_ct4 = []

model = m841


#=======================
#=====import data=======
#=======================

#=====static data=======
static_initial_cells = []
static_data_mean = []
static_data_exp1 = []
static_data_exp2 = []
static_data_exp3 = []

ct=0
with open('data/cleaned_static_data.txt', newline='', encoding='utf-8-sig') as csvfile_s:
    full_data_s = csv.reader(csvfile_s, delimiter=' ')
    for row in full_data_s:
        if ct == 0:
            static_cells = row
        elif ct > 0:
            static_initial_cells.append(float(row[0]))
            static_data_mean.append(float(row[1]))
            static_data_exp1.append(float(row[2]))
            static_data_exp2.append(float(row[3]))
            static_data_exp3.append(float(row[4]))
        ct = ct+1

#print(static_initial_cells)
#print(static_data_mean)
#print(static_data_exp1)
#print(static_data_exp2)
#print(static_data_exp3)
#exit()

#====dynamic data====
tdata = []
exp_data_DMSO = [[], [], [], []]
exp_data_10uM = [[], [], [], []]
exp_data_2_5uM = [[], [], [], []]
exp_data_625nM = [[], [], [], []]
exp_data_156_25nM = [[], [], [], []]
exp_data_39nM = [[], [], [], []]
exp_data_9_76nM = [[], [], [], []]
exp_data_2_44nM = [[], [], [], []]
exp_data_0_61nM = [[], [], [], []]
exp_data_0_15nM = [[], [], [], []]

ct=0
with open('data/cleaned_dynamic_data.txt', newline='') as csvfile_d:
    full_data_d = csv.reader(csvfile_d, delimiter=' ')
    for row in full_data_d:
        #print(row)
        if ct > 0:
            tdata.append(float(row[0]))
            exp_data_DMSO[0].append(float(row[1]))
            exp_data_DMSO[1].append(float(row[2]))
            exp_data_DMSO[2].append(float(row[3]))
            exp_data_DMSO[3].append(float(row[4]))
            exp_data_10uM[0].append(float(row[5]))
            exp_data_10uM[1].append(float(row[6]))
            exp_data_10uM[2].append(float(row[7]))
            exp_data_10uM[3].append(float(row[8]))
            exp_data_2_5uM[0].append(float(row[9]))
            exp_data_2_5uM[1].append(float(row[10]))
            exp_data_2_5uM[2].append(float(row[11]))
            exp_data_2_5uM[3].append(float(row[12]))
            exp_data_625nM[0].append(float(row[13]))
            exp_data_625nM[1].append(float(row[14]))
            exp_data_625nM[2].append(float(row[15]))
            exp_data_625nM[3].append(float(row[16]))
            exp_data_156_25nM[0].append(float(row[17]))
            exp_data_156_25nM[1].append(float(row[18]))
            exp_data_156_25nM[2].append(float(row[19]))
            exp_data_156_25nM[3].append(float(row[20]))
            exp_data_39nM[0].append(float(row[21]))
            exp_data_39nM[1].append(float(row[22]))
            exp_data_39nM[2].append(float(row[23]))
            exp_data_39nM[3].append(float(row[24]))
            exp_data_9_76nM[0].append(float(row[25]))
            exp_data_9_76nM[1].append(float(row[26]))
            exp_data_9_76nM[2].append(float(row[27]))
            exp_data_9_76nM[3].append(float(row[28]))
            exp_data_2_44nM[0].append(float(row[29]))
            exp_data_2_44nM[1].append(float(row[30]))
            exp_data_2_44nM[2].append(float(row[31]))
            exp_data_2_44nM[3].append(float(row[32]))
            exp_data_0_61nM[0].append(float(row[33]))
            exp_data_0_61nM[1].append(float(row[34]))
            exp_data_0_61nM[2].append(float(row[35]))
            exp_data_0_61nM[3].append(float(row[36]))
            exp_data_0_15nM[0].append(float(row[37]))
            exp_data_0_15nM[1].append(float(row[38]))
            exp_data_0_15nM[2].append(float(row[39]))
            exp_data_0_15nM[3].append(float(row[40]))
        ct = ct+1

#print(tdata)
#exit()
plt.plot(tdata,exp_data_DMSO[0], label = 'DMSO')
plt.plot(tdata,exp_data_0_15nM[0], label = '0.15nM')
plt.plot(tdata,exp_data_0_61nM[0], label = '0.61nM')
plt.plot(tdata,exp_data_2_44nM[0], label = '2.44nM')
plt.plot(tdata,exp_data_9_76nM[0], label = '9.76nM')
plt.plot(tdata,exp_data_39nM[0], label = '39nM')
plt.plot(tdata,exp_data_156_25nM[0], label = '156.25nM')
plt.plot(tdata,exp_data_625nM[0], label = '625nM')
plt.plot(tdata,exp_data_2_5uM[0], label = '2.5uM')
plt.plot(tdata,exp_data_10uM[0], label = '10uM')
plt.legend()
plt.title('experimental luminescence data, mean')
plt.show()
#exit()

rate_mask = np.array([True, False, False, True, True, False, True, True, True, True, False, False, False, False, False])

param_values = np.array([p.value for p in model.parameters])

solver = ScipyOdeSimulator(model, integrator='lsoda', use_analytic_jacobian=True, compiler='python')
spec_list = [str(x) for x in model.species]

init_cell = spec_list.index("Cells(NADPHoxidoreductase1=None, drug='undrugged')")
drug_idx = spec_list.index('Drug()')

lum_sim_static = []
DMSO_sim = []
lum_10um_sim = []
lum_2_5uM_sim = []
lum_625nM_sim = []
lum_156_25nM_sim = []
lum_39nM_sim = []
lum_9_76nM_sim = []
lum_2_44nM_sim = []
lum_0_61nM_sim = []
lum_0_15nM_sim = []

lum_ct_params = []
cell_string = 'H841'

cell_ct_DMSO = []
cell_ct_10uM = []
cell_ct_2_5uM = []
cell_ct_625nM = []
cell_ct_156_25nM = []
cell_ct_39nM = []
cell_ct_9_76nM = []
cell_ct_2_44nM = []
cell_ct_0_61nM = []
cell_ct_0_15nM = []
cells_drugged_DMSO = []
cells_drugged_10uM = []
cells_drugged_2_5uM = []
cells_drugged_625nM = []
cells_drugged_156_25nM = []
cells_drugged_39nM = []
cells_drugged_9_76nM = []
cells_drugged_2_44nM = []
cells_drugged_0_61nM = []
cells_drugged_0_15nM = []
cells_undrugged_DMSO = []
cells_undrugged_10uM = []
cells_undrugged_2_5uM = []
cells_undrugged_625nM = []
cells_undrugged_156_25nM = []
cells_undrugged_39nM = []
cells_undrugged_9_76nM = []
cells_undrugged_2_44nM = []
cells_undrugged_0_61nM = []
cells_undrugged_0_15nM = []



def sim_one_drug(init_drug, traj_init, new_p):
    param_values[rate_mask] = new_p
    initials = traj_init
    initials[drug_idx] = init_drug

    traj = solver.run(param_values=param_values, initials=initials, tspan=tdata)
    sim_lum = traj.expressions['Luminescence']
    sim_cell_ct = traj.observables['Live_Cells']
    sim_cells_d = traj.observables['Drugged_Live_Cells']
    sim_cells_nd = traj.observables['Undrugged_Live_Cells']
    return sim_lum, sim_cell_ct, sim_cells_d, sim_cells_nd

def sim_over_multiple_parameters(dir_string):
    with open(dir_string, newline='') as myfile:
        for line in myfile:
            new_p = [float(x) for x in line.split(',')]
            #print(new_p)
            param_values[rate_mask] = new_p

            traj_init = solver.run(param_values=param_values,
                                   initials={model.species[init_cell]: 300, model.species[drug_idx]: 0},
                                   tspan=[0, 0.5])
            initials = traj_init.species[-1]

            lum, cellct, drugged_cells, undrugged_cells = sim_one_drug(0, initials, new_p)
            DMSO_sim.append(lum)
            cell_ct_DMSO.append(cellct)
            cells_drugged_DMSO.append(drugged_cells)
            cells_undrugged_DMSO.append(undrugged_cells)

            lum, cellct, drugged_cells, undrugged_cells = sim_one_drug(10.0e-6*N_A*80e-6, initials, new_p)
            lum_10um_sim.append(lum)
            cell_ct_10uM.append(cellct)
            cells_drugged_10uM.append(drugged_cells)
            cells_undrugged_10uM.append(undrugged_cells)

            lum, cellct, drugged_cells, undrugged_cells = sim_one_drug(2.5e-6*N_A*80e-6, initials, new_p)
            lum_2_5uM_sim.append(lum)
            cell_ct_2_5uM.append(cellct)
            cells_drugged_2_5uM.append(drugged_cells)
            cells_undrugged_2_5uM.append(undrugged_cells)

            lum, cellct, drugged_cells, undrugged_cells = sim_one_drug(625e-9*N_A*80e-6, initials, new_p)
            lum_625nM_sim.append(lum)
            cell_ct_625nM.append(cellct)
            cells_drugged_625nM.append(drugged_cells)
            cells_undrugged_625nM.append(undrugged_cells)

            lum, cellct, drugged_cells, undrugged_cells = sim_one_drug(156.25e-9*N_A*80e-6, initials, new_p)
            lum_156_25nM_sim.append(lum)
            cell_ct_156_25nM.append(cellct)
            cells_drugged_156_25nM.append(drugged_cells)
            cells_undrugged_156_25nM.append(undrugged_cells)

            lum, cellct, drugged_cells, undrugged_cells = sim_one_drug(39e-9*N_A*80e-6, initials, new_p)
            lum_39nM_sim.append(lum)
            cell_ct_39nM.append(cellct)
            cells_drugged_39nM.append(drugged_cells)
            cells_undrugged_39nM.append(undrugged_cells)

            lum, cellct, drugged_cells, undrugged_cells = sim_one_drug(9.76e-9*N_A*80e-6, initials, new_p)
            lum_9_76nM_sim.append(lum)
            cell_ct_9_76nM.append(cellct)
            cells_drugged_9_76nM.append(drugged_cells)
            cells_undrugged_9_76nM.append(undrugged_cells)

            lum, cellct, drugged_cells, undrugged_cells = sim_one_drug(2.44e-9*N_A*80e-6, initials, new_p)
            lum_2_44nM_sim.append(lum)
            cell_ct_2_44nM.append(cellct)
            cells_drugged_2_44nM.append(drugged_cells)
            cells_undrugged_2_44nM.append(undrugged_cells)

            lum, cellct, drugged_cells, undrugged_cells = sim_one_drug(0.61e-9*N_A*80e-6, initials, new_p)
            lum_0_61nM_sim.append(lum)
            cell_ct_0_61nM.append(cellct)
            cells_drugged_0_61nM.append(drugged_cells)
            cells_undrugged_0_61nM.append(undrugged_cells)

            lum, cellct, drugged_cells, undrugged_cells = sim_one_drug(0.15e-9*N_A*80e-6, initials, new_p)
            lum_0_15nM_sim.append(lum)
            cell_ct_0_15nM.append(cellct)
            cells_drugged_0_15nM.append(drugged_cells)
            cells_undrugged_0_15nM.append(undrugged_cells)

    #print(sim_lum)
    return

sim_over_multiple_parameters('second_multiobj_try/gathered_params_H841.txt')


        #DMSO_sim.append(np.array(y_DMSO.all['Live_Cells']/y_DMSO.all['Live_Cells'][0]))



# plot for static data
for i in range(len(lum_sim_static)):
    if i == 0:
        plt.plot(static_initial_cells,lum_sim_static[i], 'x', label = 'simulated data')
    else:
        plt.plot(static_initial_cells, lum_sim_static[i], 'x')
plt.plot(static_initial_cells, static_data_exp1, 'o', markerfacecolor="None", label = 'experiment replicate 1', color = 'red')
plt.plot(static_initial_cells, static_data_exp2, 'o', markerfacecolor="None", label = 'experiment replicate 2', color = 'blue')
plt.plot(static_initial_cells, static_data_exp3, 'o', markerfacecolor="None", label = 'experiment replicate 3', color = 'limegreen')
plt.plot(static_initial_cells, static_data_mean, 'o', color = 'black', label = 'mean experimental data')
plt.title('static parameter estimate')
plt.xlabel('initial cell number')
plt.legend()
plt.show()
#exit()

plt.plot(tdata, cell_ct_DMSO[-1], label = 'DMSO')
plt.plot(tdata, cell_ct_0_15nM[-1], label = '0.15nM')
plt.plot(tdata, cell_ct_0_61nM[-1], label = '0.61nM')
plt.plot(tdata, cell_ct_2_44nM[-1], label = '2.44nM')
plt.plot(tdata, cell_ct_9_76nM[-1], label = '9.76nM')
plt.plot(tdata, cell_ct_39nM[-1], label = '39nM')
plt.plot(tdata, cell_ct_156_25nM[-1], label = '156.25nM')
plt.plot(tdata, cell_ct_625nM[-1], label = '625nM')
plt.plot(tdata, cell_ct_2_5uM[-1], label = '2.5uM')
plt.plot(tdata, cell_ct_10uM[-1], label = '10uM')
plt.legend()
plt.title('compare sim cell counts')
plt.show()

#print(np.log2(cell_ct_DMSO[-1]))
#exit()
plt.plot(tdata, np.log2(cell_ct_DMSO[-1]), label = 'DMSO')
plt.plot(tdata, np.log2(cell_ct_0_15nM[-1]), label = '0.15nM')
plt.plot(tdata, np.log2(cell_ct_0_61nM[-1]), label = '0.61nM')
plt.plot(tdata, np.log2(cell_ct_2_44nM[-1]), label = '2.44nM')
plt.plot(tdata, np.log2(cell_ct_9_76nM[-1]), label = '9.76nM')
plt.plot(tdata, np.log2(cell_ct_39nM[-1]), label = '39nM')
plt.plot(tdata, np.log2(cell_ct_156_25nM[-1]), label = '156.25nM')
#plt.plot(tdata, np.log2(cell_ct_625nM[-1]), label = '625nM')
plt.plot(tdata, np.log2(cell_ct_2_5uM[-1]), label = '2.5uM')
plt.plot(tdata, np.log2(cell_ct_10uM[-1]), label = '10uM')
plt.legend()
plt.ylabel('log2')
plt.title('compare sim cell counts')
plt.show()



for i in range(len(cell_ct_DMSO)):
    #plt.plot(tdata,cell_ct[i], 'x', label = 'simulated data')
    if i == 0:
        plt.plot(tdata, cells_undrugged_DMSO[i], 'x', label='undrugged')
        plt.plot(tdata, cells_drugged_DMSO[i], 'o', label = 'drugged')
        plt.plot(tdata, cell_ct_DMSO[i], label = 'total')
    else:
        plt.plot(tdata, cells_undrugged_DMSO[i], 'x')
        plt.plot(tdata, cells_drugged_DMSO[i], 'o')
        plt.plot(tdata, cell_ct_DMSO[i])
plt.title('DMSO - cell counts')
plt.legend()
plt.show()

#plot for DMSO dynamic data
for i in range(len(DMSO_sim)):
    if i == 0:
        #plt.plot(tdata, DMSO_sim[i], 'x', label = 'simulated data')
        plt.plot(tdata, DMSO_sim[i], label = 'simulated data')
    else:
        #plt.plot(tdata, DMSO_sim[i], 'x')
        plt.plot(tdata, DMSO_sim[i])
#plt.plot(tdata, exp_data_DMSO[1], 'o', markerfacecolor="None", label = 'experiment replicate 1', color = 'red')
#plt.plot(tdata, exp_data_DMSO[2], 'o', markerfacecolor="None", label = 'experiment replicate 2', color = 'blue')
#plt.plot(tdata, exp_data_DMSO[3], 'o', markerfacecolor="None", label = 'experiment replicate 3', color = 'limegreen')
plt.plot(tdata, exp_data_DMSO[0], 'o', color = 'black', label = 'mean experimental data')
plt.legend()
plt.xlabel('time')
plt.ylabel('luminescence')
plt.title('dynamic parameter estimate: DMSO')
plt.show()
#exit()



for i in range(len(cell_ct_10uM)):
    #plt.plot(tdata,cell_ct[i], 'x', label = 'simulated data')
    if i == 0:
        plt.plot(tdata, cells_undrugged_10uM[i], 'x', label='undrugged')
        plt.plot(tdata, cells_drugged_10uM[i], 'o', label = 'drugged')
        plt.plot(tdata, cell_ct_10uM[i], label = 'total')
    else:
        plt.plot(tdata, cells_undrugged_10uM[i], 'x')
        plt.plot(tdata, cells_drugged_10uM[i], 'o')
        plt.plot(tdata, cell_ct_10uM[i])
plt.title('10uM - cell counts')
plt.legend()
plt.show()

#plot for 10uM dynamic data
for i in range(len(lum_10um_sim)):
    if i == 0:
        #plt.plot(tdata, lum_10um_sim[i], 'x', label = 'simulated data')
        plt.plot(tdata, lum_10um_sim[i], label = 'simulated data')
    else:
        #plt.plot(tdata, lum_10um_sim[i], 'x')
        plt.plot(tdata, lum_10um_sim[i])
#plt.plot(tdata, exp_data_10uM[1], 'o', markerfacecolor="None", label = 'experiment replicate 1', color = 'red')
#plt.plot(tdata, exp_data_10uM[2], 'o', markerfacecolor="None", label = 'experiment replicate 2', color = 'blue')
#plt.plot(tdata, exp_data_10uM[3], 'o', markerfacecolor="None", label = 'experiment replicate 3', color = 'limegreen')
plt.plot(tdata, exp_data_10uM[0], 'o', color = 'black', label = 'mean experimental data')
plt.legend()
plt.xlabel('time')
plt.ylabel('luminescence')
plt.title('dynamic parameter estimate: 10uM')
plt.show()
#exit()



for i in range(len(cell_ct_2_5uM)):
    #plt.plot(tdata,cell_ct[i], 'x', label = 'simulated data')
    if i == 0:
        plt.plot(tdata, cells_undrugged_2_5uM[i], 'x', label='undrugged')
        plt.plot(tdata, cells_drugged_2_5uM[i], 'o', label = 'drugged')
        plt.plot(tdata, cell_ct_2_5uM[i], label = 'total')
    else:
        plt.plot(tdata, cells_undrugged_2_5uM[i], 'x')
        plt.plot(tdata, cells_drugged_2_5uM[i], 'o')
        plt.plot(tdata, cell_ct_2_5uM[i])
plt.title('2.5uM - cell counts')
plt.legend()
plt.show()


#plot for 2.5uM dynamic data
for i in range(len(lum_2_5uM_sim)):
    if i == 0:
        #plt.plot(tdata, lum_2_5uM_sim[i], 'x', label = 'simulated data')
        plt.plot(tdata, lum_2_5uM_sim[i], label = 'simulated data')
    else:
        #plt.plot(tdata, lum_2_5uM_sim[i], 'x')
        plt.plot(tdata, lum_2_5uM_sim[i])
#plt.plot(tdata, exp_data_2_5uM[1], 'o', markerfacecolor="None", label = 'experiment replicate 1', color = 'red')
#plt.plot(tdata, exp_data_2_5uM[2], 'o', markerfacecolor="None", label = 'experiment replicate 2', color = 'blue')
#plt.plot(tdata, exp_data_2_5uM[3], 'o', markerfacecolor="None", label = 'experiment replicate 3', color = 'limegreen')
plt.plot(tdata, exp_data_2_5uM[0], 'o', color = 'black', label = 'mean experimental data')
plt.legend()
plt.xlabel('time')
plt.ylabel('luminescence')
plt.title('dynamic parameter estimate: 2.5uM')
plt.show()
#exit()




for i in range(len(cell_ct_625nM)):
    #plt.plot(tdata,cell_ct[i], 'x', label = 'simulated data')
    if i == 0:
        plt.plot(tdata, cells_undrugged_625nM[i], 'x', label='undrugged')
        plt.plot(tdata, cells_drugged_625nM[i], 'o', label = 'drugged')
        plt.plot(tdata, cell_ct_625nM[i], label = 'total')
    else:
        plt.plot(tdata, cells_undrugged_625nM[i], 'x')
        plt.plot(tdata, cells_drugged_625nM[i], 'o')
        plt.plot(tdata, cell_ct_625nM[i])
plt.title('625nM - cell counts')
plt.legend()
plt.show()

#plot for 625nM dynamic data
for i in range(len(lum_625nM_sim)):
    if i == 0:
        #plt.plot(tdata, lum_625nM_sim[i], 'x', label = 'simulated data')
        plt.plot(tdata, lum_625nM_sim[i], label = 'simulated data')
    else:
        #plt.plot(tdata, lum_625nM_sim[i], 'x')
        plt.plot(tdata, lum_625nM_sim[i])
#plt.plot(tdata, exp_data_625nM[1], 'o', markerfacecolor="None", label = 'experiment replicate 1', color = 'red')
#plt.plot(tdata, exp_data_625nM[2], 'o', markerfacecolor="None", label = 'experiment replicate 2', color = 'blue')
#plt.plot(tdata, exp_data_625nM[3], 'o', markerfacecolor="None", label = 'experiment replicate 3', color = 'limegreen')
plt.plot(tdata, exp_data_625nM[0], 'o', color = 'black', label = 'mean experimental data')
plt.legend()
plt.xlabel('time')
plt.ylabel('luminescence')
plt.title('dynamic parameter estimate: 625nM')
plt.show()
#exit()



for i in range(len(cell_ct_156_25nM)):
    #plt.plot(tdata,cell_ct[i], 'x', label = 'simulated data')
    if i == 0:
        plt.plot(tdata, cells_undrugged_156_25nM[i], 'x', label='undrugged')
        plt.plot(tdata, cells_drugged_156_25nM[i], 'o', label = 'drugged')
        plt.plot(tdata, cell_ct_156_25nM[i], label = 'total')
    else:
        plt.plot(tdata, cells_undrugged_156_25nM[i], 'x')
        plt.plot(tdata, cells_drugged_156_25nM[i], 'o')
        plt.plot(tdata, cell_ct_156_25nM[i])
plt.title('156.25nM - cell counts')
plt.legend()
plt.show()


#plot for 156.25nM dynamic data
for i in range(len(lum_156_25nM_sim)):
    if i == 0:
        #plt.plot(tdata, lum_156_25nM_sim[i], 'x', label = 'simulated data')
        plt.plot(tdata, lum_156_25nM_sim[i], label = 'simulated data')
    else:
        #plt.plot(tdata, lum_156_25nM_sim[i], 'x')
        plt.plot(tdata, lum_156_25nM_sim[i])
#plt.plot(tdata, exp_data_156_25nM[1], 'o', markerfacecolor="None", label = 'experiment replicate 1', color = 'red')
#plt.plot(tdata, exp_data_156_25nM[2], 'o', markerfacecolor="None", label = 'experiment replicate 2', color = 'blue')
#plt.plot(tdata, exp_data_156_25nM[3], 'o', markerfacecolor="None", label = 'experiment replicate 3', color = 'limegreen')
plt.plot(tdata, exp_data_156_25nM[0], 'o', color = 'black', label = 'mean experimental data')
plt.legend()
plt.xlabel('time')
plt.ylabel('luminescence')
plt.title('dynamic parameter estimate: 156.25nM')
plt.show()
#exit()



for i in range(len(cell_ct_39nM)):
    #plt.plot(tdata,cell_ct[i], 'x', label = 'simulated data')
    if i == 0:
        plt.plot(tdata, cells_undrugged_39nM[i], 'x', label='undrugged')
        #plt.plot(tdata, cells_drugged_39nM[i], 'o', label = 'drugged')
        plt.plot(tdata, cell_ct_39nM[i], label = 'total')
    else:
        plt.plot(tdata, cells_undrugged_39nM[i], 'x')
        #plt.plot(tdata, cells_drugged_39nM[i], 'o')
        plt.plot(tdata, cell_ct_39nM[i])
plt.title('39nM - cell counts')
plt.legend()
plt.show()


#plot for 39nM dynamic data
for i in range(len(lum_39nM_sim)):
    if i == 0:
        #plt.plot(tdata, lum_39nM_sim[i], 'x', label = 'simulated data')
        plt.plot(tdata, lum_39nM_sim[i], label = 'simulated data')
    else:
        #plt.plot(tdata, lum_39nM_sim[i], 'x')
        plt.plot(tdata, lum_39nM_sim[i])
#plt.plot(tdata, exp_data_39nM[1], 'o', markerfacecolor="None", label = 'experiment replicate 1', color = 'red')
#plt.plot(tdata, exp_data_39nM[2], 'o', markerfacecolor="None", label = 'experiment replicate 2', color = 'blue')
#plt.plot(tdata, exp_data_39nM[3], 'o', markerfacecolor="None", label = 'experiment replicate 3', color = 'limegreen')
plt.plot(tdata, exp_data_39nM[0], 'o', color = 'black', label = 'mean experimental data')
plt.legend()
plt.xlabel('time')
plt.ylabel('luminescence')
plt.title('dynamic parameter estimate: 39nM')
plt.show()
print(lum_39nM_sim[1][0])
#exit()



for i in range(len(cell_ct_9_76nM)):
    #plt.plot(tdata,cell_ct[i], 'x', label = 'simulated data')
    if i == 0:
        plt.plot(tdata, cells_undrugged_9_76nM[i], 'x', label='undrugged')
        #plt.plot(tdata, cells_drugged_9_76nM[i], 'o', label = 'drugged')
        plt.plot(tdata, cell_ct_9_76nM[i], label = 'total')
    else:
        plt.plot(tdata, cells_undrugged_9_76nM[i], 'x')
        #splt.plot(tdata, cells_drugged_9_76nM[i], 'o')
        plt.plot(tdata, cell_ct_9_76nM[i])
plt.title('9.76nM - cell counts')
plt.legend()
plt.show()


#plot for 9.76nM dynamic data
for i in range(len(lum_9_76nM_sim)):
    if i == 0:
        #plt.plot(tdata, lum_9_76nM_sim[i], 'x', label = 'simulated data')
        plt.plot(tdata, lum_9_76nM_sim[i], label = 'simulated data')
    else:
        #plt.plot(tdata, lum_9_76nM_sim[i], 'x')
        plt.plot(tdata, lum_9_76nM_sim[i])
#plt.plot(tdata, exp_data_9_76nM[1], 'o', markerfacecolor="None", label = 'experiment replicate 1', color = 'red')
#plt.plot(tdata, exp_data_9_76nM[2], 'o', markerfacecolor="None", label = 'experiment replicate 2', color = 'blue')
#plt.plot(tdata, exp_data_9_76nM[3], 'o', markerfacecolor="None", label = 'experiment replicate 3', color = 'limegreen')
plt.plot(tdata, exp_data_9_76nM[0], 'o', color = 'black', label = 'mean experimental data')
plt.legend()
plt.xlabel('time')
plt.ylabel('luminescence')
plt.title('dynamic parameter estimate: 9.76nM')
plt.show()
print(lum_9_76nM_sim[0][0])
#exit()


for i in range(len(cell_ct_2_44nM)):
    #plt.plot(tdata,cell_ct[i], 'x', label = 'simulated data')
    if i == 0:
        plt.plot(tdata, cells_undrugged_2_44nM[i], 'x', label='undrugged')
        plt.plot(tdata, cells_drugged_2_44nM[i], 'o', label = 'drugged')
        plt.plot(tdata, cell_ct_2_44nM[i], label = 'total')
    else:
        plt.plot(tdata, cells_undrugged_2_44nM[i], 'x')
        plt.plot(tdata, cells_drugged_2_44nM[i], 'o')
        plt.plot(tdata, cell_ct_2_44nM[i])
plt.title('2.44nM - cell counts')
plt.legend()
plt.show()


#plot for 2.44nM dynamic data
for i in range(len(lum_2_44nM_sim)):
    if i == 0:
        #plt.plot(tdata, lum_2_44nM_sim[i], 'x', label = 'simulated data')
        plt.plot(tdata, lum_2_44nM_sim[i], label = 'simulated data')
    else:
        #plt.plot(tdata, lum_2_44nM_sim[i], 'x')
        plt.plot(tdata, lum_2_44nM_sim[i])
#plt.plot(tdata, exp_data_2_44nM[1], 'o', markerfacecolor="None", label = 'experiment replicate 1', color = 'red')
#plt.plot(tdata, exp_data_2_44nM[2], 'o', markerfacecolor="None", label = 'experiment replicate 2', color = 'blue')
#plt.plot(tdata, exp_data_2_44nM[3], 'o', markerfacecolor="None", label = 'experiment replicate 3', color = 'limegreen')
plt.plot(tdata, exp_data_2_44nM[0], 'o', color = 'black', label = 'mean experimental data')
plt.legend()
plt.xlabel('time')
plt.ylabel('luminescence')
plt.title('dynamic parameter estimate: 2.44nM')
plt.show()
print(lum_2_44nM_sim[0][0])
#exit()



for i in range(len(cell_ct_0_61nM)):
    #plt.plot(tdata,cell_ct[i], 'x', label = 'simulated data')
    if i == 0:
        plt.plot(tdata, cells_undrugged_0_61nM[i], 'x', label='undrugged')
        plt.plot(tdata, cells_drugged_0_61nM[i], 'o', label = 'drugged')
        plt.plot(tdata, cell_ct_0_61nM[i], label = 'total')
    else:
        plt.plot(tdata, cells_undrugged_0_61nM[i], 'x')
        plt.plot(tdata, cells_drugged_0_61nM[i], 'o')
        plt.plot(tdata, cell_ct_0_61nM[i])
plt.title('0.61nM - cell counts')
plt.legend()
plt.show()


#plot for 0.61nM dynamic data
for i in range(len(lum_0_61nM_sim)):
    if i == 0:
        #plt.plot(tdata, lum_0_61nM_sim[i], 'x', label = 'simulated data')
        plt.plot(tdata, lum_0_61nM_sim[i], label = 'simulated data')
    else:
        #plt.plot(tdata, lum_0_61nM_sim[i], 'x')
        plt.plot(tdata, lum_0_61nM_sim[i])
#plt.plot(tdata, exp_data_0_61nM[1], 'o', markerfacecolor="None", label = 'experiment replicate 1', color = 'red')
#plt.plot(tdata, exp_data_0_61nM[2], 'o', markerfacecolor="None", label = 'experiment replicate 2', color = 'blue')
#plt.plot(tdata, exp_data_0_61nM[3], 'o', markerfacecolor="None", label = 'experiment replicate 3', color = 'limegreen')
plt.plot(tdata, exp_data_0_61nM[0], 'o', color = 'black', label = 'mean experimental data')
plt.legend()
plt.xlabel('time')
plt.ylabel('luminescence')
plt.title('dynamic parameter estimate: 0.61nM')
plt.show()
print(lum_0_61nM_sim[0][0])
#exit()




for i in range(len(cell_ct_0_15nM)):
    #plt.plot(tdata,cell_ct[i], 'x', label = 'simulated data')
    if i == 0:
        plt.plot(tdata, cells_undrugged_0_15nM[i], 'x', label='undrugged')
        plt.plot(tdata, cells_drugged_0_15nM[i], 'o', label = 'drugged')
        plt.plot(tdata, cell_ct_0_15nM[i], label = 'total')
    else:
        plt.plot(tdata, cells_undrugged_0_15nM[i], 'x')
        plt.plot(tdata, cells_drugged_0_15nM[i], 'o')
        plt.plot(tdata, cell_ct_0_15nM[i])
plt.title('0.15nM - cell counts')
plt.legend()
plt.show()

#plot for 0.15nM dynamic data
for i in range(len(lum_0_15nM_sim)):
    if i == 0:
        #plt.plot(tdata, lum_0_15nM_sim[i], 'x', label = 'simulated data')
        plt.plot(tdata, lum_0_15nM_sim[i], label = 'simulated data')
    else:
        #plt.plot(tdata, lum_0_15nM_sim[i], 'x')
        plt.plot(tdata, lum_0_15nM_sim[i])
#plt.plot(tdata, exp_data_0_15nM[1], 'o', markerfacecolor="None", label = 'experiment replicate 1', color = 'red')
#plt.plot(tdata, exp_data_0_15nM[2], 'o', markerfacecolor="None", label = 'experiment replicate 2', color = 'blue')
#plt.plot(tdata, exp_data_0_15nM[3], 'o', markerfacecolor="None", label = 'experiment replicate 3', color = 'limegreen')
plt.plot(tdata, exp_data_0_15nM[0], 'o', color = 'black', label = 'mean experimental data')
plt.legend()
plt.xlabel('time')
plt.ylabel('luminescence')
plt.title('dynamic parameter estimate: 0.15nM')
plt.show()
print(lum_0_15nM_sim[0][0])
#exit()


plt.plot(tdata, exp_data_DMSO[0], 'o', color = 'red')
plt.plot(tdata, DMSO_sim[-1],  color = 'red', label = 'DMSO')
plt.plot(tdata, exp_data_0_15nM[0], 'o', color = 'blue')
plt.plot(tdata, lum_0_15nM_sim[-1],  color = 'blue', label = '0.15nM')
plt.plot(tdata, exp_data_0_61nM[0], 'o', color = 'green')
plt.plot(tdata, lum_0_61nM_sim[-1],  color = 'green', label = '0.61nM')
plt.plot(tdata, exp_data_2_44nM[0], 'o', color = 'gold')
plt.plot(tdata, lum_2_44nM_sim[-1],  color = 'gold', label = '2.44nM')
plt.plot(tdata, exp_data_9_76nM[0], 'o', color = 'purple')
plt.plot(tdata, lum_9_76nM_sim[-1],  color = 'purple', label = '9.76nM')
plt.plot(tdata, exp_data_39nM[0], 'o', color = 'turquoise')
plt.plot(tdata, lum_39nM_sim[-1],  color = 'turquoise', label = '39nM')
plt.plot(tdata, exp_data_156_25nM[0], 'o', color = 'hotpink')
plt.plot(tdata, lum_156_25nM_sim[-1],  color = 'hotpink', label = '156.25nM')
plt.plot(tdata, exp_data_625nM[0], 'o', color = 'dodgerblue')
plt.plot(tdata, lum_625nM_sim[-1],  color = 'dodgerblue', label = '625nM')
plt.plot(tdata, exp_data_2_5uM[0], 'o', color = 'darkorange')
plt.plot(tdata, lum_2_5uM_sim[-1],  color = 'darkorange', label = '2.5uM')
plt.plot(tdata, exp_data_10uM[0], 'o', color = 'lawngreen')
plt.plot(tdata, lum_10um_sim[-1],  color = 'lawngreen', label = '10uM')
plt.title('exp data vs sim lum')
plt.legend()
plt.show()


for i in range(len(lum_0_15nM_sim)):
    if i == 0:
        #plt.plot(tdata, lum_0_15nM_sim[i], 'x', label = 'simulated data')
        plt.plot(tdata, lum_0_15nM_sim[i], label = '0.15nM')
        plt.plot(tdata, lum_10um_sim[i], 'x', label = '10uM')
    else:
        #plt.plot(tdata, lum_0_15nM_sim[i], 'x')
        plt.plot(tdata, lum_0_15nM_sim[i])
        plt.plot(tdata, lum_10um_sim[i], 'x')
plt.show()