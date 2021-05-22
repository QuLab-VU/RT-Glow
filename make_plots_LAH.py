from model_H841_k3d import model
# from model_H841 import model
from pysb.simulator import ScipyOdeSimulator
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import N_A
from get_data_LAH import get_data
import re
import os
from string import ascii_uppercase, ascii_letters

# experimental data
FILENAME = 'data/H1048_DMS53_Combined_Lum_and_cell_counts.csv' #'data/H841_Lum_only.csv'
CELL_LINE = 'NCIH1048' #'H841' 
DRUG_NAME = 'TAK-901'

exp_data = get_data(FILENAME, CELL_LINE, DRUG_NAME)
drug_conc_all = [0, 10000, 2500, 625, 156, 39.1, 9.77, 2.44, 0.61, 0.15]

'''
# plots
for i in range(len(drug_conc_all)):
    plt.figure('cells (%g)' % drug_conc_all[i])
    #     plt.plot([tdata]*len(cell_count_reps), cell_count_reps, 'o')
    plt.plot(tdata[i], cell_count_mean[i], 'ko')
    plt.xlabel('time (h)')
    plt.ylabel('cell count')
    plt.figure('lum (%g)' % drug_conc_all[i])
    #     plt.plot([tdata]*len(lum_rlu_reps), lum_rlu_reps, 's')
    plt.plot(tdata[i], lum_rlu_mean[i], 'ks')
    plt.xlabel('time (h)')
    plt.ylabel('luminescence (RLU)')
'''
# exp_data = np.genfromtxt('data/cleaned_dynamic_data.txt', names=True)
# 
# col_names_all = ['DMSO_mean', '10uM_mean', 'cleaned_2.5uM_mean', 'cleaned_625nM_mean', 
#                  'cleaned_156.25nM_mean', 'cleaned_39nM_mean', 'cleaned_9.76nM_mean', 
#                  'cleaned_2.44nM_mean', 'cleaned_0.61nM_mean', 'cleaned_0.15nM_mean']
# 
# drug_conc_all = [0, 10e-6, 2.5e-6, 625e-9, 156.25e-9, 39e-9, 9.76e-9, 2.44e-9, 0.61e-9, 0.15e-9]

#####
drug_conc = [[conc*1e-9] for conc in drug_conc_all] 
col_names = [[str(conc)+'nM'] for conc in drug_conc_all] 
'''
for i in range(len(col_names)):
    for conc,col in zip(drug_conc[i],col_names[i]):
        plt.figure('lum ({:.2f} nM)'.format(conc*1e9))
        plt.plot(exp_data['tdata_%s' % col], exp_data['lum_%s_mean' % col], 
                 'o', color='k', label='exp data')
        plt.xlabel('time (h)')
        plt.ylabel('luminescence (RLU)')
        plt.legend(loc=0,ncol=2)
plt.show()
quit()
'''
#####    

dir = '%s_%s' % (CELL_LINE, DRUG_NAME) #'H841_TAK901' #'H1048_TAK901'
n_iter = 100
suffix = '_w_k3d_iter_%d' % n_iter
max_num_psets = 9 #20 # default should be # of PSO particles (20 in this case)

outdir = '%s_%s' % (CELL_LINE, DRUG_NAME)    
OUTFILE = open(os.path.join(outdir, 
                            '%s_iter_%d_cell_counts.csv' 
                            % (model.name, n_iter)),"w")
OUTFILE.write('upid,well,cell.count,time,cell.line,drug1,drug1.conc,drug1.units\n')

# plate_row = 0 # A
plate_col = 1
for i in range(len(drug_conc)):
    
    # remove dots ('.') from data column names 
#     for j in range(len(col_names[i])):
#         col_names[i][j] = col_names[i][j].replace('.', '')

    # fitted parameter values
    infile = os.path.join(dir, 'run%d%s.csv' % (i,suffix))
    fit_params = np.genfromtxt(infile, dtype=None, delimiter=',', names=True)
    # print(fit_params.dtype.names)

    psets = [fp for fp in fit_params if fp['iter'] == (n_iter-1) and fp['fitness'] < 1]
    if len(psets) > max_num_psets:
        psets = psets[:max_num_psets]
    cmap = plt.cm.jet
    colors = [cmap(f) for f in np.linspace(0,1,len(psets))]
    
    # initial cell count
    cell_init = 300
    
    # parameter names and values
    pnames = np.array([p.name for p in model.parameters])
    indices = [np.where(pnames == name)[0][0] for name in fit_params.dtype.names 
               if len(np.where(pnames == name)[0]) > 0]
    pvals = np.array([p.value for p in model.parameters])
    
    sim = ScipyOdeSimulator(model, verbose=True)
    
    sp_names = [str(sp) for sp in model.species]
    cell_idx = sp_names.index("Cells(NADPHoxidoreductase1=None, drug='undrugged')")
    drug_idx = sp_names.index("Drug()")
        
    # loop over parameter sets
    plate_row = 0 # A
    for pset,color in zip(psets, colors):
        
        # remove 1st and last 2 entries ('index', 'fitness', 'particle')    
        pvals[indices] = 10**np.array(list(pset)[1:-2])
                
        # Pre-simulation (no drug)
        def presim(param_values, cell_init, t_max, n_pts):
            
            tspan = np.linspace(0,t_max,n_pts)
            x = sim.run(tspan=tspan, param_values=param_values, 
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
            for k in range(len(slope)-1):
                if slope[k+1] > slope[k]:
                    break
            equil_idx = k
            t_equil = tspan[equil_idx]
            tspan_equil = tspan[:equil_idx+1] - t_equil
            
            # make sure equilibrium has been reached
            # if not, return a large cost function value
            equil_reached = True
            if (k == len(slope)-2):
                equil_reached = False
        
            ###### Plots ######
            '''
            for drug in drug_conc[i]:
                plt.figure('lum ({:.2f} nM)'.format(drug*1e9))
                plt.plot(0, x.expressions['Luminescence'][:equil_idx+1][-1], 
                         'x', color=color, ms=12, mew=2)
                plt.plot(tspan_equil, x.expressions['Luminescence'][:equil_idx+1], 
                         color=color, lw=2)
                
                plt.figure('cells ({:.2f} nM)'.format(drug*1e9))
    #             plt.plot(0, np.log2(x.observables['Live_Cells'][:equil_idx+1][-1]/cell_init), 
    #                      'x', color=color, ms=12, mew=2)
    #             plt.plot(tspan_equil, np.log2(x.observables['Live_Cells'][:equil_idx+1]/cell_init), 
    #                      color=color, lw=2, label='sse %.3f' % pset['fitness'])
                plt.plot(0, x.observables['Live_Cells'][:equil_idx+1][-1], 
                         'x', color=color, ms=12, mew=2)
                plt.plot(tspan_equil, x.observables['Live_Cells'][:equil_idx+1], 
                         color=color, lw=2)
            '''
            return equil_reached, x.species[equil_idx]
                
        equil_reached, initials = presim(pvals, cell_init, 1, 10001)
                
        ###### Add drug ######
        
        tspan = np.linspace(0,5,501)
        
        for conc,col in zip(drug_conc[i],col_names[i]):
        
            print(conc)
        
            initials[drug_idx] = conc*N_A*80e-6
            x = sim.run(tspan=tspan, param_values=pvals, initials=initials)
            
            # luminescence
            plt.figure('lum ({:.2f} nM)'.format(conc*1e9))
            plt.plot(tspan, x.expressions['Luminescence'], 
                     color=color, lw=2, label='%.3f' % pset['fitness'])
            plt.xlabel('time (d)')
            plt.ylabel('luminescence')
#             plt.legend(loc=0,ncol=2)
            
            # cell counts
            plt.figure('cells ({:.2f} nM)'.format(conc*1e9))
            plt.plot(tspan, x.observables['Live_Cells'], 
                     color=color, lw=2, label='fit %.3f' % pset['fitness'])
            plt.ylabel('cell count')            
#             cells_t0 = x.observables['Live_Cells'][0]
#             plt.plot(tspan, np.log2(x.observables['Live_Cells']/cells_t0), 
#                      color=color, lw=2, label='%.3f' % pset['fitness'])
#             plt.ylabel('doublings')
            plt.xlabel('time (d)')
#             plt.legend(loc=0,ncol=2)
            
            ### output to file ###
            upid = '%s_iter_%d' % (model.name, n_iter)
            well = '%c%.2d' % (ascii_uppercase[plate_row], plate_col)
            cell_count = x.observables['Live_Cells'] # array
            time = tspan*24. # array (in hours)
            cell_line = 'v'+CELL_LINE # 'v' for 'virtual'
            drug1 = DRUG_NAME
            drug1_conc = conc
            drug1_units = 'M'

            for t in range(len(time)):
#                 print(upid, well, cell_count[t], time[t]*24., cell_line, 
#                       drug1, drug1_conc, drug1_units)
                if (col != '2500nM'): # TEMPORARY
                    OUTFILE.write('%s,%s,%g,%g,%s,%s,%g,%s\n' %
                          (upid, well, cell_count[t], time[t], cell_line,
                          drug1, drug1_conc, drug1_units))
        
        # increase row number in virtual plate
        plate_row += 1                
#     # increase column number in virtual plate
#     plate_col += 1
            
#   experimental data
    for conc,col in zip(drug_conc[i],col_names[i]):
        
        title = '%s + %s: %g nM' % (CELL_LINE, DRUG_NAME, conc*1e9)
        
        # luminescence
        plt.figure('lum ({:.2f} nM)'.format(conc*1e9))
        plt.title(title)
        plt.plot(exp_data['tdata_%s' % col], exp_data['lum_%s_mean' % col], 
                 'o', color='k', label='exp data')
        plt.legend(loc=0, ncol=2, title='Lum fit (sum-of-squared errors)')
#                    bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        
        # cell counts
        plt.figure('cells ({:.2f} nM)'.format(conc*1e9))
        plt.title(title)
        plt.plot(exp_data['tdata_%s' % col], exp_data['cells_%s_mean' % col], 
                 's', color='k', label='exp data')
#         cells_t0 = exp_data['cells_%s_mean' % col][0]
#         plt.plot(exp_data['tdata_%s' % col], np.log2(exp_data['cells_%s_mean' % col]/cells_t0), 
#                  's', color='k', label='exp data')
        plt.legend(loc=0, ncol=2, title='Lum fit (sum-of-squared errors)') 
#                    bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout() #pad, h_pad, w_pad, rect)

        ### output to file ###
        upid = '%s_iter_%d' % (model.name, n_iter)
        well = '%c%.2d' % (ascii_uppercase[plate_row], plate_col)
        cell_count = exp_data['cells_%s_mean' % col] # array
        time = exp_data['tdata_%s' % col]*24. # array (in hours)
        cell_line = CELL_LINE # 'v' for 'virtual'
        drug1 = DRUG_NAME
        drug1_conc = conc
        drug1_units = 'M'
            
        for t in range(len(time)):
            OUTFILE.write('%s,%s,%g,%g,%s,%s,%g,%s\n' %
                  (upid, well, cell_count[t], time[t], cell_line,
                  drug1, drug1_conc, drug1_units))

    # increase column number in virtual plate
    plate_col += 1

OUTFILE.close()
plt.show()






