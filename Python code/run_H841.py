from model_H841 import model
import numpy as np
from pysb.simulator import ScipyOdeSimulator
from simplepso.pso import PSO
from math import sqrt
from scipy.constants import N_A
from numpy.core.numeric import inf
import sys
import csv

from scipy import stats
import matplotlib.pyplot as plt
import math as m


#=======================
#=====import data=======
#=======================

static_initial_cells= []
static_measured_lum = []

ct=0
with open('data/cleaned_static_data.txt', newline='') as csvfile_s:
    full_data_s = csv.reader(csvfile_s, delimiter=' ')
    for row in full_data_s:
        if ct>0:
            static_initial_cells.append(float(row[0]))
            static_measured_lum.append(float(row[1]))
        ct = ct+1
#print(static_initial_cells)
#print(static_measured_lum)


tdata = []
exp_DMSO = []
exp_10uM = []
exp_2_5uM =[]
exp_625nM = []
exp_156_25nM = []
exp_39nM = []
exp_9_76nM = []
exp_2_44nM = []
exp_0_61nM =[]
exp_0_15nM = []

ct=0
with open('data/cleaned_dynamic_data.txt', newline='') as csvfile_s:
    full_data_s = csv.reader(csvfile_s, delimiter=' ')
    for row in full_data_s:
        if ct>0:
            tdata.append(float(row[0]))
            exp_DMSO.append(float(row[1]))
            exp_10uM.append(float(row[5]))
            exp_2_5uM.append(float(row[9]))
            exp_625nM.append(float(row[13]))
            exp_156_25nM.append(float(row[17]))
            exp_39nM.append(float(row[21]))
            exp_9_76nM.append(float(row[25]))
            exp_2_44nM.append(float(row[29]))
            exp_0_61nM.append(float(row[33]))
            exp_0_15nM.append(float(row[37]))
        ct = ct+1



#=======================
#=====prepare PSO=======
#=======================

solver = ScipyOdeSimulator(model, integrator = 'lsoda', use_analytic_jacobian=True, compiler='python')

spec_list = [str(x) for x in model.species]
init_cell = spec_list.index("Cells(NADPHoxidoreductase1=None, drug='undrugged')")
drug_idx = spec_list.index('Drug()')

param_values = np.array([p.value for p in model.parameters])

# parameters: k3(T), k5(F), k6(F), km(T), kcat(T), cc(F), kon(T), koff(T), k5d(T), k6d(T), gain(F), cells0(F), nanoluc0(F), prosub0(F), drug0(F)
# rate mask: k3(T), km(T), kcat(T), kon(T), koff(T), k5d(T), k6d(T)
rate_mask = np.array([True, False, False, True, True, False, True, True, True, True, False, False, False, False, False])

original_values = np.array([p.value for p in model.parameters])
k5 = original_values[1]
k6 = original_values[2]

# We search in log10 space for the parameters
starting_position = np.log10(param_values[rate_mask])


#--------function for the static cell data---------------

def obj_static(position):
    param_values[rate_mask] = 10 ** position.copy()

    lum_ct = []
    tspan = np.linspace(0, 0.5, 10)
    for i in static_initial_cells:
        traj = solver.run(param_values=param_values, initials={model.species[init_cell]: i, model.species[drug_idx]: 0}, tspan=tspan)
        lum_ct.append(traj.expressions['Luminescence'][-1])

    return np.sum((np.array(static_measured_lum) - np.array(lum_ct)) ** 2) / (np.sum(np.array(static_measured_lum) ** 2))


#--------functions for the dynamic cell data---------------

def obj_dynamic(position, cellcount, drug_amount, drug_data):
    param_values[rate_mask] = 10 ** position.copy()
    #tspan_init = np.linspace(0, 0.5, 2)
    traj_init = solver.run(param_values=param_values,
                      initials={model.species[init_cell]: cellcount, model.species[drug_idx]: 0},
                      tspan=[0,0.5])
    #print(traj_init.expressions['Luminescence'])
    #exit()
    initials = traj_init.species[-1]   #start luminescence value with value getting from the initial simulation (pass in the substrate concentration)
    initials[drug_idx]=drug_amount
    traj = solver.run(param_values=param_values,
                      initials= initials,
                      tspan=tdata)
    lum_traj = traj.expressions['Luminescence']

    return np.sum((np.array(drug_data) - lum_traj) **2)/(np.sum(np.array(drug_data)**2))


# ----- PSO -----

def costfunction(parameter):
    if k5 < parameter[5]:  #k5 needs to keep larger than k5d     k3(T), km(T), kcat(T), kon(T), koff(T), k5d(T), k6d(T)
        return inf,
    if k6 > parameter[6]:  #k6 needs to keep smaller than k6d
        return inf,
    error_static = 1*obj_static(parameter)

    error0 = obj_dynamic(parameter, 300, 0, exp_DMSO)
    error1 = obj_dynamic(parameter, 300, 10.0e-6*N_A*80e-6, exp_10uM)
    error2 = obj_dynamic(parameter, 300, 2.5e-6*N_A*80e-6, exp_2_5uM)
    error3 = obj_dynamic(parameter, 300, 625e-9*N_A*80e-6, exp_625nM)
    error4 = obj_dynamic(parameter, 300, 156.25e-9*N_A*80e-6, exp_156_25nM)
    error5 = obj_dynamic(parameter, 300, 39e-9*N_A*80e-6, exp_39nM)
    error6 = obj_dynamic(parameter, 300, 9.76e-9*N_A*80e-6, exp_9_76nM)
    error7 = obj_dynamic(parameter, 300, 2.44e-9*N_A*80e-6, exp_2_44nM)
    error8 = obj_dynamic(parameter, 300, 0.61e-9*N_A*80e-6, exp_0_61nM)
    error9 = obj_dynamic(parameter, 300, 0.15e-9*N_A*80e-6, exp_0_15nM)

    total_err = error_static +error0+error1+error2+error3+error4+error5+error6+error7+error8+error9
    #total_err = error0
    if np.isnan(total_err):
        return inf,
    else:
        return total_err,


def run_pso(run, iterations, bd):
    pso = PSO(save_sampled=False, verbose = True, shrink_steps = False, num_proc=14)
    pso.set_cost_function(costfunction)
    pso.set_start_position(starting_position)
    pso.set_bounds(bd)
    pso.set_speed(-.1,.1)

    pso.run(num_particles=200, num_iterations=iterations, stop_threshold=1e-5)
    #print('best pos: ', pso.best.pos)
    print('history ', pso.history)
    print('run ', run)
    print('fit ', pso.best.fitness)
    print('all fitness ', pso.values)
    np.savetxt("H841_params_"+run+".txt", 10**pso.history, delimiter=",")
    np.savetxt("H841_fit_" + run + ".txt", pso.values, delimiter=",")


def main():
    
    arg_lst = []
    for arg in sys.argv[1:]:
        #print(arg)
        arg_lst.append(arg)

    #print(arg_lst)
    #exit()
    run = arg_lst[0]
    #data_file = arg_lst[1]
    iterations = int(arg_lst[1])
    bounds = float(arg_lst[2])

    run_pso(run = run, iterations = iterations, bd = bounds)
    '''
    run = str(1)
    iterations = 100
    bounds = 4
    run_pso(run = run, iterations = iterations, bd = bounds)
    '''

if '__main__' == __name__:
    main()








