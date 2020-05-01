from model_H841 import model
import numpy as np
from pysb.simulator import ScipyOdeSimulator
from simplepso.pso import PSO
from math import sqrt
from scipy.constants import N_A
from numpy.core.numeric import inf
import sys
import csv
import statistics

from scipy import stats
import matplotlib.pyplot as plt
import math as m

#-------import data-------

# --- DMSO data ---

tdata_DMSO = []
exp_data_DMSO = []


ct=0
with open('mKate2counts_TotHour.csv', newline='') as csvfile:
    full_data = csv.reader(csvfile, delimiter=',')
    for row in full_data:
        if ct > 0:
            #print(row)
            tdata_DMSO.append(float(row[19]))
            #exp_data.append(float(row[20]))  DMSO Data for H1048
            #exp_data.append(float(row[21]))  DMSO Data for H1048
            exp_data_DMSO.append(float(row[22]))  #DMSO Data for H841
        ct = ct+1

tdata_DMSO = [x/24 for x in tdata_DMSO]
print(tdata_DMSO)

# --- experimental data ---

tdata_d = []
exp_data_10uM = []
exp_data_2_5uM = []
exp_data_625nM = []
exp_data_156_25nM = []
exp_data_39nM = []
exp_data_9_76nM = []
exp_data_2_44nM = []
exp_data_0_61nM = []
exp_data_0_15nM = []

ct=0
with open('H841_TAK901_Lum.csv', newline='') as csvfile_d:
    full_data_d = csv.reader(csvfile_d, delimiter=',')
    for row in full_data_d:
        #print(row)
        if ct > 0:
            tdata_d.append(float(row[0]))
            exp_data_10uM.append(float(row[1]))
            exp_data_2_5uM.append(float(row[2]))
            exp_data_625nM.append(float(row[3]))
            exp_data_156_25nM.append(float(row[4]))
            exp_data_39nM.append(float(row[5]))
            exp_data_9_76nM.append(float(row[6]))
            exp_data_2_44nM.append(float(row[7]))
            exp_data_0_61nM.append(float(row[8]))
            exp_data_0_15nM.append(float(row[9]))
        ct = ct+1

tdata_drug = [x/24 for x in tdata_d]

# static data

static_cells5000 = []
static_cells2500 = []
static_cells1250 = []
static_cells625 = []
static_cells312 = []
static_cells156 = []
static_cells78 = []
static_cells39 = []
static_cells0 = []

ct=0
with open('H841_Static_Lum_5k.csv', newline='', encoding='utf-8-sig') as csvfile_s:
    full_data_s = csv.reader(csvfile_s, delimiter=',')
    for row in full_data_s:
        if ct == 0:
            static_cells = row
        elif ct > 0:
            static_cells5000.append(float(row[0]))
            static_cells2500.append(float(row[1]))
            static_cells1250.append(float(row[2]))
            static_cells625.append(float(row[3]))
            static_cells312.append(float(row[4]))
            static_cells156.append(float(row[5]))
            static_cells78.append(float(row[6]))
            static_cells39.append(float(row[7]))
            static_cells0.append(float(row[8]))
        ct = ct+1
static_initial_ct = [float(x) for x in static_cells]
print(static_cells)
#exit()

static_cells_mean = [statistics.mean(static_cells0), statistics.mean(static_cells39), statistics.mean(static_cells78),
                     statistics.mean(static_cells156), statistics.mean(static_cells312), statistics.mean(static_cells625),
                     statistics.mean(static_cells1250), statistics.mean(static_cells2500), statistics.mean(static_cells5000)]

static_log_cells = [m.log2(x) for x in static_cells_mean]
#plt.plot(static_log_cells)
#plt.show()

slope_static, intercept_static, r_static, p_static, std_static = stats.linregress(range(len(static_log_cells)),static_log_cells)

#plt.plot(range(len(static_log_cells)), static_log_cells, 'x')
#plt.plot(range(len(static_log_cells)), intercept_static+range(len(static_log_cells))*slope_static)
#plt.show()
#exit()

#print(slope_static, r_static)
#exit()

#print(tdata_DMSO)
#print(tdata_drug)

#--------prepare PSO stuff-------------------
#solver = ScipyOdeSimulator(model, tspan=tdata_drug, integrator='lsoda',use_analytic_jacobian=True, compiler='python', cleanup=False, verbose=1)
#solver = ScipyOdeSimulator(model, tspan=tdata_drug, integrator='vode', compiler='python', cleanup=False, verbose=1)
#y = solver.run()
#print(y.expressions['Luminescence'])
#exit()
#solver_DMSO = ScipyOdeSimulator(model, tspan=tdata_DMSO, integrator='vode',use_analytic_jacobian=True, compiler='python')
#solver_drug = ScipyOdeSimulator(model, tspan=tdata_DMSO, integrator='vode',use_analytic_jacobian=True, compiler='python')
#solver_drug = ScipyOdeSimulator(model, tspan=tdata_drug, integrator='lsoda', integrator_options={atol:1e-6, rtol:1e-6},use_analytic_jacobian=True, compiler='python')


# --- PSO preparations ---

solver = ScipyOdeSimulator(model, integrator = 'lsoda', use_analytic_jacobian=True, compiler='python')
#spec_list = model.species
spec_list = [str(x) for x in model.species]
init_cell = spec_list.index("Cells(NADPHoxidoreductase1=None, drug='undrugged')")
drug_idx = spec_list.index('Drug()')
#exit()
param_values = np.array([p.value for p in model.parameters])

# parameters: k3(T), k5(F), k6(F), km(T), kcat(T), cc(F), kon(T), koff(T), k5d(T), k6d(T), gain(F), cells0(F), nanoluc0(F), prosub0(F), drug0(F)
# rate mask: k3(T), km(T), kcat(T), kon(T), koff(T), k5d(T), k6d(T)
rate_mask = np.array([True, False, False, True, True, False, True, True, True, True, False, False, False, False, False])

original_values = np.array([p.value for p in model.parameters])
k5 = original_values[1]
k6 = original_values[2]

# We search in log10 space for the parameters
starting_position = np.log10(param_values[rate_mask])


#--------function for the dynamic cell data---------------

def obj_static(position):
    param_values[rate_mask] = 10 ** position.copy()

    lum_ct = []
    rev = [x for x in reversed(static_initial_ct)]
    #print(rev)
    #exit()
    tspan = np.linspace(0, 0.5, 100)
    for i in rev:
        #print("i: ",i)
        traj = solver.run(param_values=param_values, initials={model.species[init_cell]: i}, tspan=tspan)
        #plt.figure()
        #plt.plot(tspan, traj.expressions['Luminescence'])
        lum_ct.append(traj.expressions['Luminescence'][-1])
    #plt.plot(rev,lum_ct, 'x')
    #plt.show()
    #exit()
    slope_sim, intercept_sim, r_sim, p_sim, std_sim = stats.linregress(range(len(lum_ct)),lum_ct)

    return (slope_static - slope_sim)**2
    #err = (slope_static - slope_sim)**2 #+ (r_static - r_sim)**2
    #return sqrt(err)


#--------functions for the dynamic cell data---------------

def obj_DMSO(position):
    param_values[rate_mask] = 10 ** position.copy()
    #param_values[14] = 0
    #traj = solver_DMSO.run(param_values=param_values)
    #print(model.species)
    #exit()
    traj = solver.run(param_values=param_values, initials={model.species[drug_idx]: 0}, tspan=tdata_DMSO)
    lum_traj = traj.expressions['Luminescence']

    return np.sum((exp_data_DMSO - lum_traj) **2)
    #err = np.sum((exp_data_DMSO - lum_traj) **2)
    #return sqrt(err)

def obj_10uM(position):
    param_values[rate_mask] = 10 ** position.copy()
    #param_values[14] = 10.0e-6*N_A*80e-6
    traj = solver.run(param_values=param_values, initials={model.species[drug_idx]: 10.0e-6*N_A*80e-6}, tspan=tdata_drug)
    lum_traj = traj.expressions['Luminescence']

    return np.sum((exp_data_10uM - lum_traj) **2)
    #err = np.sum((exp_data_10uM - lum_traj) **2)
    #return sqrt(err)

def obj_2_5uM(position):
    param_values[rate_mask] = 10 ** position.copy()
    #param_values[14] = 2.5e-6*N_A*80e-6
    traj = solver.run(param_values=param_values, initials={model.species[drug_idx]: 2.5e-6*N_A*80e-6}, tspan=tdata_drug)
    lum_traj = traj.expressions['Luminescence']

    return np.sum((exp_data_2_5uM - lum_traj) **2)
    #err = np.sum((exp_data_2_5uM - lum_traj) **2)
    #return sqrt(err)

def obj_625nM(position):
    param_values[rate_mask] = 10 ** position.copy()
    #param_values[14] = 625e-9*N_A*80e-6
    traj = solver.run(param_values=param_values, initials={model.species[drug_idx]: 625e-9*N_A*80e-6}, tspan=tdata_drug)
    lum_traj = traj.expressions['Luminescence']

    return np.sum((exp_data_625nM - lum_traj) **2)
    #err = np.sum((exp_data_625nM - lum_traj) **2)
    #return sqrt(err)

def obj_156_25nM(position):
    param_values[rate_mask] = 10 ** position.copy()
    #param_values[14] = 156.25e-9*N_A*80e-6
    traj = solver.run(param_values=param_values, initials={model.species[drug_idx]: 156.25e-9*N_A*80e-6}, tspan=tdata_drug)
    lum_traj = traj.expressions['Luminescence']

    return np.sum((exp_data_156_25nM - lum_traj) **2)
    #err = np.sum((exp_data_156_25nM - lum_traj) **2)
    #return sqrt(err)

def obj_39nM(position):
    param_values[rate_mask] = 10 ** position.copy()
    #param_values[14] = 39e-9*N_A*80e-6
    traj = solver.run(param_values=param_values, initials={model.species[drug_idx]: 39e-9*N_A*80e-6}, tspan=tdata_drug)
    lum_traj = traj.expressions['Luminescence']

    return np.sum((exp_data_39nM - lum_traj) **2)
    #err = np.sum((exp_data_39nM - lum_traj) **2)
    #return sqrt(err)

def obj_9_76nM(position):
    param_values[rate_mask] = 10 ** position.copy()
    #param_values[14] = 9.76e-9*N_A*80e-6
    traj = solver.run(param_values=param_values, initials={model.species[drug_idx]: 9.76e-9*N_A*80e-6}, tspan=tdata_drug)
    lum_traj = traj.expressions['Luminescence']

    return np.sum((exp_data_9_76nM - lum_traj) **2)
    #err = np.sum((exp_data_9_76nM - lum_traj) **2)
    #return sqrt(err)

def obj_2_44nM(position):
    param_values[rate_mask] = 10 ** position.copy()
    #param_values[14] = 2.44e-9*N_A*80e-6
    traj = solver.run(param_values=param_values, initials={model.species[drug_idx]: 2.44e-9*N_A*80e-6}, tspan=tdata_drug)
    lum_traj = traj.expressions['Luminescence']

    return np.sum((exp_data_2_44nM - lum_traj) **2)
    #err = np.sum((exp_data_2_44nM - lum_traj) **2)
    #return sqrt(err)

def obj_0_61nM(position):
    param_values[rate_mask] = 10 ** position.copy()
    #param_values[14] = 0.61e-9*N_A*80e-6
    traj = solver.run(param_values=param_values, initials={model.species[drug_idx]: 0.61e-9*N_A*80e-6}, tspan=tdata_drug)
    lum_traj = traj.expressions['Luminescence']

    return np.sum((exp_data_0_61nM - lum_traj) **2)
    #err = np.sum((exp_data_0_61nM - lum_traj) **2)
    #return sqrt(err)

def obj_0_15nM(position):
    param_values[rate_mask] = 10 ** position.copy()
    #param_values[14] = 0.15e-9*N_A*80e-6
    traj = solver.run(param_values=param_values, initials={model.species[drug_idx]: 0.15e-9*N_A*80e-6}, tspan=tdata_drug)
    lum_traj = traj.expressions['Luminescence']

    return np.sum((exp_data_0_15nM - lum_traj) **2)
    #err = np.sum((exp_data_0_15nM - lum_traj) **2)
    #return sqrt(err)

# ----- PSO -----

def costfunction(parameter):
    if k5 > parameter[5]:  #k5 > k5d     k3(T), km(T), kcat(T), kon(T), koff(T), k5d(T), k6d(T)
        return inf,
    if k6 < parameter[6]:  #k6 < k6d
        return inf,
    error_static = obj_static(parameter)
    error0 = obj_DMSO(parameter)
    error1 = obj_10uM(parameter)
    error2 = obj_2_5uM(parameter)
    error3 = obj_625nM(parameter)
    error4 = obj_156_25nM(parameter)
    error5 = obj_39nM(parameter)
    error6 = obj_9_76nM(parameter)
    error7 = obj_2_44nM(parameter)
    error8 = obj_0_61nM(parameter)
    error9 = obj_0_15nM(parameter) 

    total_err = error_static+error0+error1+error2+error3+error4+error5+error6+error7+error8+error9
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

    pso.run(num_particles=100, num_iterations=iterations, stop_threshold=1e-5)
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
    bounds = int(arg_lst[2])

    run_pso(run = run, iterations = iterations, bd = bounds)
    '''
    run = str(1)
    iterations = 100
    bounds = 4
    run_pso(run = run, iterations = iterations, bd = bounds)
    '''

if '__main__' == __name__:
    main()








