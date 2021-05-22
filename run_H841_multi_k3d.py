# from model_H841_k3d import model
from model_H841 import model
from simplepso.pso import PSO
from model_H841_k3d import model
import os
import pandas as pd
#os.environ['BNGPATH'] = r'C:\Users\pinojc\Desktop\BioNetGen-2.5.1'
import numpy as np
import matplotlib.pyplot as plt
from pysb.simulator import ScipyOdeSimulator #, BngSimulator
from scipy.constants import N_A
import sys
import csv
import re
import os
from get_data_LAH import get_data

#=======================
#=====import data=======
#=======================
'''
static_initial_cells = []
static_measured_lum = []
ct = 0
with open('data/cleaned_static_data.txt', newline='') as csvfile_s:
    full_data_s = csv.reader(csvfile_s, delimiter=' ')
    for row in full_data_s:
        if ct > 0:
            static_initial_cells.append(float(row[0]))
            static_measured_lum.append(float(row[1]))
        ct = ct+1
'''
# exp_data = pd.read_csv('data/cleaned_dynamic_data.txt', delimiter='\s')


'''
def get_drug_conc_and_column_names():
    # get drug concentrations from column names
    conc = []
    names = []
    for col in exp_data.columns:
        match = re.search("DMSO_", col)
        if match:
            conc.append(0)
            names.append(col)
        else:
            match = re.search("(\d+\.*\d*)(\w+)_", col)
            if match:
                units = match.group(2)
                if units == 'uM':
                    scale = 1e-6
                elif units == 'nM':
                    scale = 1e-9
                else:
                    print('Error: Concentration units \'%s\' not recognized.' % units)
                    quit()
                conc.append(scale*float(match.group(1)))
                names.append(col)
    # print(drug_conc)
    # quit()
    return conc, names       
'''
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
# generate_equations(model)
spec_list = [str(sp) for sp in model.species]
cell_idx = spec_list.index("Cells(NADPHoxidoreductase1=None, drug='undrugged')")
drug_idx = spec_list.index('Drug()')

param_names = np.array([p.name for p in model.parameters])
param_values = np.array([p.value for p in model.parameters])

k5_val = param_values[list(param_names).index('K5')]
k6_val = param_values[list(param_names).index('K6')]

# parameters: k3(T), k5(F), k6(F), km(T), kcat(T), cc(F), kon(T), koff(T), k3d(T), k5d(T), k6d(T), gain(F), cells0(F), nanoluc0(F), prosub0(F), drug0(F)
# rate mask: k3(T), km(T), kcat(T), kon(T), koff(T), k3d(T), k5d(T), k6d(T)
free_params = ['K3', 'Km', 'Kcat', 'kon', 'koff', 'K3d', 'K5d', 'K6d']
rate_mask = np.array([p.name in free_params for p in model.parameters])
# rate_mask_old = np.array([True, False, False, True, True, False, True, True, True, True, True, False, False, False, False, False])

k5d_idx = list(param_names[rate_mask]).index('K5d')
k6d_idx = list(param_names[rate_mask]).index('K6d')

# We search in log10 space for the parameters
starting_params = np.log10(param_values[rate_mask])

def presim(param_values, cell_init, t_max, n_pts):

    tspan = np.linspace(0,t_max,n_pts)
        
    x = solver.run(tspan=tspan, param_values=param_values, 
        initials={
            model.species[cell_idx]: cell_init,
            model.species[drug_idx]: 0}
        )
    
    # find equilibrium time point
    lum = x.expressions['Luminescence']
    delta_t = tspan[1]-tspan[0]
    slope = (np.roll(lum,-1)[:-1]-lum[:-1])/delta_t
    # slope starts out large, declines, and then increases again
    # define equilibrium as first point at which slope begins increasing
    for i in range(len(slope)-1):
        if slope[i+1] > slope[i]:
            break
    equil_idx = i
    t_equil = tspan[equil_idx]
    tspan_equil = tspan[:equil_idx+1] - t_equil
    
    # make sure equilibrium has been reached
    # if not, return a large cost function value
    equil_reached = True
    if (i == len(slope)-2):
        equil_reached = False

    return equil_reached, x.species[equil_idx]

# --------function for the static cell data---------------

def obj_static(position):
    
    param_values[rate_mask] = 10**position.copy()

    lum_ct = []
    tspan = np.linspace(0, 0.5, 10)
    for i in static_initial_cells:
        traj = solver.run(param_values=param_values,
                          initials={model.species[cell_idx]: i,
                                    model.species[drug_idx]: 0}, 
                          tspan=tspan)
        lum_ct.append(traj.expressions['Luminescence'][-1])

    return np.sum((np.array(static_measured_lum) - np.array(lum_ct)) ** 2) / (
        np.sum(np.array(static_measured_lum) ** 2))

# --------functions for the dynamic cell data---------------

def obj_dynamic(position, cellcount, drug_amount, t_data, drug_data, tol):

    verbose = False
    if verbose:
        print('compute initials for drug {} using a lum start of {}'
              ''.format(drug_amount, drug_data[0]))

    param_values[rate_mask] = 10**position.copy()
    
    # run a pre-simulation w/ no drug
    equil_reached, initials = presim(param_values, cellcount, 1, 10001)
    if not equil_reached:
        return 1e8
    
    # now add the drug and run another simulation
    initials[drug_idx] = drug_amount
    traj = solver.run(
        param_values=param_values,
        initials=initials,
        tspan=t_data #exp_data['tdata'].values
    )
    lum = traj.expressions['Luminescence']
    sse = np.sum(((drug_data - lum)/drug_data) ** 2)
    
    return sse

# ----- PSO -----

def costfunction(parameter):
    
    # k5 must be larger than k5d (division rate)
    if k5_val < 10**parameter[k5d_idx]:
        return 1e8
    # k6 must be smaller than k6d (death rate)
    if k6_val > 10**parameter[k6d_idx]:
        return 1e8
    
    tol = 1e-5
    
    # if 'drug_conc' is not a list, make it a list
    try:
        iter(drug_conc)
    except TypeError:
        conc = [drug_conc]
    else:
        conc = drug_conc
    
    # if 'data_column' is not a list, make it a list
    try:
        iter(data_column)
    except TypeError:
        col = [data_column]
    else:
        col = data_column
        
    # error check
    if len(conc) != len(col):
        print('Error: length of \'drug_conc\' (%d) and \'data_column\' (%d)' + 
              'must be equal' % (len(conc), len(col)))
        quit()
    
    err = [None]*len(conc)
    for j in range(len(conc)):
        err[j] = obj_dynamic(parameter, 300, conc[j]*N_A*80e-6, 
                             exp_data['tdata_%s' % col[j]],
                             exp_data['lum_%s_mean' % col[j]], tol)

#     error0 = obj_dynamic(parameter, 300, 0.15e-9*N_A*80e-6, 
#                          exp_data['cleaned_0.15nM_mean'], tol)
    
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
#     total_err = error0
    total_err = sum(err)
    
    if np.isnan(total_err):
        return 1e8
    else:
        return total_err

def convert_to_flat_array(optimizer, model):
    # convert to array
    history = np.array(optimizer.all_history)
    fitness = np.array(optimizer.all_fitness)
#     print(fitness)
#     print(fitness.shape)
#     print(history.shape)
#     quit()
    col_names = list(param_names[rate_mask]) + ['fitness']
#     print(col_names)

    # convert to pandas dataframe
    all_df = []
    for i in range(history.shape[1]):
        pos = history[:, i]
        fit_value = fitness[:, i]
#         print(fit_value)
        stacked = np.column_stack([pos, fit_value])
        new_df = pd.DataFrame(stacked, columns=col_names)
        new_df['particle'] = i
        all_df.append(new_df)
    
    return pd.concat(all_df)

def run_pso(run, iterations, bd, outdir = '', suffix=''):
    
    pso = PSO(save_sampled=True, verbose=True, shrink_steps=False)
    
    #pso.set_cost_function(costfunction)
    pso.set_start_position(starting_params)
    pso.set_bounds(bd)
    pso.set_speed(-0.1, 0.1)

    pso.run(num_particles=20, num_iterations=iterations, stop_threshold=1e-5,
            cost_function=costfunction, max_iter_no_improv=500,
            num_processors=20, save_samples=True)
    
    param_sets = convert_to_flat_array(pso, model)
    
    #print(param_sets)
    outfile = os.path.join(outdir, 'run%d%s.csv' % (run,suffix))
    param_sets.to_csv(outfile, index_label='iter')
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
        iterations = int(arg_lst[1])
        bounds = float(arg_lst[2])
    else:
        run = 0
        iterations = 100
        bounds = 1
    
    FILENAME = 'data/H1048_DMS53_Combined_Lum_and_cell_counts.csv' #'data/H841_Lum_only.csv'
    CELL_LINE = 'NCIH1048' #'H841'
    DRUG_NAME = 'TAK-901'
    
    global exp_data 
    exp_data = get_data(FILENAME, CELL_LINE, DRUG_NAME)
    
#     drug_conc_all, col_names_all = get_drug_conc_and_column_names()
    drug_conc_all = [0, 10000, 2500, 625, 156, 39.1, 9.77, 2.44, 0.61, 0.15] # nM

#     drug_conc_samples = [[drug] for drug in drug_conc_all]    
#     col_names_samples = [[col] for col in col_names_all] 

    drug_conc_samples = [[conc*1e-9] for conc in drug_conc_all] 
    col_names_samples = [[str(conc)+'nM'] for conc in drug_conc_all] 
    
    # make output directory
    outdir = '%s_%s' % (CELL_LINE, DRUG_NAME) #'H841_TAK901' #'H1048_TAK901'
    os.makedirs(outdir, exist_ok=True)
    
    global drug_conc
    global data_column
    for i in range(len(drug_conc_samples)):
        drug_conc = drug_conc_samples[i] 
        data_column = col_names_samples[i] 
        print(drug_conc, data_column)
        run_pso(run=i, iterations=iterations, bd=bounds, 
                outdir = outdir, suffix='_no_k3d_iter_%d' % iterations)

if '__main__' == __name__:
    main()








