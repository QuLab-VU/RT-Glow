from simplepso.pso import PSO
from model_H841 import model
import os
import pandas as pd
#os.environ['BNGPATH'] = r'C:\Users\pinojc\Desktop\BioNetGen-2.5.1'
import numpy as np
from pysb.simulator import ScipyOdeSimulator, BngSimulator
from math import sqrt
from scipy.constants import N_A
from numpy.core.numeric import inf
from pysb.bng import generate_equations
import sys
import csv
from pysb.pathfinder import set_path
#from pysb.generator.bng import BngGenerator

from scipy import stats
import matplotlib.pyplot as plt
import math as m

#=======================
#=====import data=======
#=======================

static_initial_cells= []
static_measured_lum = []

ct = 0
with open('data/cleaned_static_data.txt', newline='') as csvfile_s:
    full_data_s = csv.reader(csvfile_s, delimiter=' ')
    for row in full_data_s:
        if ct > 0:
            static_initial_cells.append(float(row[0]))
            static_measured_lum.append(float(row[1]))
        ct = ct+1
#print(static_initial_cells)
#print(static_measured_lum)

exp_data = pd.read_csv('data/cleaned_dynamic_data.txt', delimiter='\s')

#=======================
#=====prepare PSO=======
#=======================

solver = ScipyOdeSimulator(
    model,
    integrator='lsoda',
    use_analytic_jacobian=True,
    compiler='cython',
    integrator_options={'atol': 1e-6, 'rtol': 1e-6}
)

# solver = BngSimulator(model, verbose=False, method='ode')
generate_equations(model)
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


def presim(param_values, cellcount, drug_seed, endpoints, refinements):
    # check, if parameter set could not find a good fit for the initial concentration:

    percent_range = drug_seed * 0.1

    #def check_index(stop, spacing):
    tspan = np.linspace(0, endpoints, int(refinements))
    traj_init = solver.run(

        param_values=param_values,
        initials={
            model.species[init_cell]: cellcount,
            model.species[drug_idx]: 0
        },
        tspan=tspan,
    )
    lum = traj_init.expressions['Luminescence']



    # find equilibrium
    diff = (np.roll(lum,-1)-lum)/(tspan[1]-tspan[0])
    diff = abs(diff/abs(lum))

    diff_loc = diff < 2
    diff_loc_1 = np.where(diff_loc)[0]

    '''
    print('diff value: ', diff[diff_loc][0], diff_loc_1[0])

    plt.plot(lum, 'x')
    plt.title('luminescence')
    plt.show()
    plt.plot(diff, 'x')
    plt.title('diff')
    plt.show()
    '''

    #lum_locate = abs(lum - drug_seed) < percent_range
    #locations = lum[lum_locate]

    #print()
    if len(lum[diff_loc]) == 0:
        return [], np.inf
    else:
        return traj_init.species[diff_loc][0], abs(lum[diff_loc][0] - drug_seed)/drug_seed


'''
        if len(locations) > 0:
            return traj_init.species[lum_locate][0]
        else:
            return None
        

    # value = check_index(5, 1e7)
    # if value is not None:
    #     return value
    for ep in endpoints:
        for r in refinements:
            value = check_index(ep, r)
            if value is not None:
                return value
    return []
'''


# --------function for the static cell data---------------

def obj_static(position):
    param_values[rate_mask] = 10 ** position.copy()

    lum_ct = []
    tspan = np.linspace(0, 0.5, 10)
    for i in static_initial_cells:
        traj = solver.run(param_values=param_values,
                          initials={model.species[init_cell]: i,
                                    model.species[drug_idx]: 0}, tspan=tspan)
        lum_ct.append(traj.expressions['Luminescence'][-1])

    return np.sum((np.array(static_measured_lum) - np.array(lum_ct)) ** 2) / (
        np.sum(np.array(static_measured_lum) ** 2))


# --------functions for the dynamic cell data---------------

def obj_dynamic(position, cellcount, drug_amount, drug_data, tol):
    param_values[rate_mask] = 10 ** position.copy()
    verbose = False
    if verbose:
        print('compute initials for drug {} using a lum start of {}'
              ''.format(drug_amount, drug_data[0]))

    initials, fit = presim(
        param_values, cellcount, drug_data[0],
        1, 1e4
        #[1, 0.5, 0.1, 0.01, 0.001, 0.0001, 1e-5, 1e-6, 1e-7, 1e-8],
        #[1e4, 1e5, 1e6]
    )


    if len(initials) == 0:
        if verbose:
            print("NOOOOOOOOO")
        return 1e8

    initials[drug_idx] = drug_amount

    traj = solver.run(
        param_values=param_values,
        initials=initials,
        tspan=exp_data['tdata'].values
    )
    lum_traj = traj.expressions['Luminescence']
    error = np.sum((drug_data - lum_traj) ** 2) / np.sum(drug_data ** 2)
    return error


# ----- PSO -----

def costfunction(parameter):
    tol = 1e-5
    # k5 needs to keep larger than k5d k3(T), km(T), kcat(T), kon(T), koff(T),
    # k5d(T), k6d(T)
    if k5 < parameter[5]:
        return 1e8
    # k6 needs to keep smaller than k6d
    if k6 > parameter[6]:
        return 1e8
    # error_static = 1*obj_static(parameter)
    # print('error static: ', error_static)

    error0 = obj_dynamic(parameter, 300, 0, exp_data['DMSO_mean'], tol)

    # 'DMSO_mean', '10uM_mean', 'cleaned_2.5uM_mean',
    # 'cleaned_625nM_mean', 'cleaned_156.25nM_mean', 'cleaned_39nM_mean',
    # 'cleaned_9.76nM_mean', 'cleaned_2.44nM_mean', 'cleaned_0.61nM_mean',
    # 'cleaned_0.15nM_mean'
    # print('param values: ', parameter)
    # print(exp_10uM[0])

    # replace all the explicit lists with the values as below

    # error1 = obj_dynamic(parameter, 300, 10.0e-6*N_A*80e-6,
    #                      exp_data['10uM_mean'], tol)

    # print('param values: ', parameter)
    # print(exp_2_5uM[0])
    # error2 = obj_dynamic(parameter, 300, 2.5e-6*N_A*80e-6, exp_2_5uM, exp_2_5uM[0], tol)
    # exit()
    # print('param values: ', parameter)
    # error3 = obj_dynamic(parameter, 300, 625e-9*N_A*80e-6, exp_625nM, exp_625nM[0], tol)
    # print('error 3: ', error3)
    # error4 = obj_dynamic(parameter, 300, 156.25e-9*N_A*80e-6, exp_156_25nM, exp_156_25nM[0], tol)
    # print('error 4: ', error4)
    # error5 = obj_dynamic(parameter, 300, 39e-9*N_A*80e-6, exp_39nM, exp_39nM[0], tol)
    # print('error 5: ', error5)
    # error6 = obj_dynamic(parameter, 300, 9.76e-9*N_A*80e-6, exp_9_76nM, exp_9_76nM[0], tol)
    # print('error 6: ', error6)
    # error7 = obj_dynamic(parameter, 300, 2.44e-9*N_A*80e-6, exp_2_44nM, exp_2_44nM[0], tol)
    # print('error 7: ', error7)
    # error8 = obj_dynamic(parameter, 300, 0.61e-9*N_A*80e-6, exp_0_61nM, exp_0_61nM[0], tol)
    # print('error 8: ', error8)
    # error9 = obj_dynamic(parameter, 300, 0.15e-9*N_A*80e-6, exp_0_15nM, exp_0_15nM[0], tol)
    # print('error 9: ', error9)
    # exit()

    # total_err = error_static +error0+error1+error2+error3+error4+error5+error6+error7+error8+error9
    # total_err = error3+error4+error5+error6+error7+error8+error9
    total_err = error0
    if np.isnan(total_err):
        return 1e8
    else:
        return total_err

def convert_to_flat_array(optimizer, model):
    # convert to array
    history = np.array(optimizer.all_history)
    fitness = np.array(optimizer.all_fitness)
    # create name masks
    rate_params = model.parameters_rules()
    rate_mask = np.array([p in rate_params for p in model.parameters])
    param_values = np.array([p.value for p in model.parameters])
    param_names = np.array([p.name for p in model.parameters])
    col_names = list(param_names[rate_mask]) + ['fitness']
    # convert to pandas dataframe
    all_df = []
    for i in range(history.shape[1]):
        pos = history[:, 0]
        fit_value = fitness[:, 0]
        stacked = np.column_stack([pos, fit_value])
        new_df = pd.DataFrame(stacked, columns=col_names)
        new_df['particle'] = i
        all_df.append(new_df)
    return pd.concat(all_df)


def run_pso(run, iterations, bd):
    pso = PSO(save_sampled=False, verbose=True, shrink_steps=False)
    #pso.set_cost_function(costfunction)
    pso.set_start_position(starting_position)
    pso.set_bounds(bd)
    pso.set_speed(-.1, .1)

    pso.run(num_particles=20, num_iterations=iterations, stop_threshold=1e-5,
            cost_function=costfunction, max_iter_no_improv=50,
            num_processors=20, save_samples=True)

    param_sets = convert_to_flat_array(pso, model)
    #print(param_sets)
    param_sets.to_csv('run'+run+'.csv')
    # print('best pos: ', pso.best.pos)
    #hist_flat = [10**item for sublist in pso.all_history for item in sublist]
    #hist_flat = pso.all_history.ravel()
    #hist_flat = pso.all_history.flatten()
    #print('history ', hist_flat)
    #print('run ', run)
    #print('fit: ', pso.all_fitness)
    #fit_flat = [item for sublist in pso.all_fitness for item in sublist]
    #print('fit ', fit_flat)
    #print('all fitness ', pso.values)
    #exit()
    #np.savetxt("H841_params_" + str(run) + ".txt", hist_flat, delimiter=",")
    #np.savetxt("H841_fit_" + str(run) + ".txt", fit_flat, delimiter=",")


def main():
    
    arg_lst = []
    for arg in sys.argv[1:]:
        arg_lst.append(arg)

    if len(arg_lst) == 3:
        run = arg_lst[0]
        data_file = arg_lst[1]
        iterations = int(arg_lst[1])
        bounds = float(arg_lst[2])
    else:
        run = 0
        iterations = 100
        bounds = 1
    run_pso(run=run, iterations=iterations, bd=bounds)
    '''
    run = str(1)
    iterations = 100
    bounds = 4
    run_pso(run = run, iterations = iterations, bd = bounds)
    '''

if '__main__' == __name__:
    main()








