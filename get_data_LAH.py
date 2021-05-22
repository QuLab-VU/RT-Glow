import numpy as np
import pandas as pd

def get_data(infile, cell_line, drug_name):

    raw_data = np.genfromtxt(infile, names=True, dtype=None, delimiter=',', encoding='UTF-8')
#     print(raw_data)
#     print(raw_data.dtype.names)

    exp_data = np.array([d for d in raw_data 
                         if d['Cell_Line'] == cell_line 
                         and d['Drug'] == drug_name])
    
    wells = np.unique([d['Well'] for d in exp_data])
    well_data = []
    for well in wells:
        well_data.append([])
        well_data[-1] = np.array([d for d in exp_data if d['Well'] == well])
    well_data = np.array(well_data)
     
    drug_conc = np.unique([d['Drug_Conc'] for d in exp_data])
    replicates = {}
    for conc in drug_conc:
        replicates[conc] = []
        for well in well_data:
            if well[0]['Drug_Conc'] == conc:
                replicates[conc].append(well)
        replicates[conc] = np.array(replicates[conc])
    
    tdata = []
    cell_count_reps = []
    cell_count_mean = []
    lum_rlu_reps = []
    lum_rlu_mean = []
    for key in replicates.keys():
        cell_count_reps.append([])
        lum_rlu_reps.append([])
        for rep in replicates[key]:
            cell_count_reps[-1].append(rep['Cell_Count'])
            lum_rlu_reps[-1].append(rep['RLU'])
        tdata.append(replicates[key][0]['TotHour']/24.) # time in days
        cell_count_mean.append(np.mean(cell_count_reps[-1], axis=0))    
        lum_rlu_mean.append(np.mean(lum_rlu_reps[-1], axis=0))
    
    data_dict = {}
    for i in range(len(drug_conc)):
        # time points
        data_dict['tdata_%gnM' % (drug_conc[i]*1e9)] = tdata[i]
        # cell counts (mean + replicates)
        data_dict['cells_%gnM_mean' % (drug_conc[i]*1e9)] = cell_count_mean[i]
        for j in range(len(cell_count_reps[i])):
            data_dict['cells_%gnM_rep%d' % (drug_conc[i]*1e9, j)] = cell_count_reps[i][j]
        # luminescence (mean + replicates)
        data_dict['lum_%gnM_mean' % (drug_conc[i]*1e9)] = lum_rlu_mean[i]
        for j in range(len(lum_rlu_reps[i])):
            data_dict['lum_%gnM_rep%d' % (drug_conc[i]*1e9, j)] = lum_rlu_reps[i][j]
    
    return pd.DataFrame(data_dict)

# data = get_data('data/H841_Lum_only.csv', 'H841', 'TAK-901')
# 
# print(data)
# print(data.columns)
# # for col in data.columns:
# #     print(data[col])
# #     print()
# 
# drug_conc_all = [0, 10000, 2500, 625, 156, 39.1, 9.77, 2.44, 0.61, 0.15]
