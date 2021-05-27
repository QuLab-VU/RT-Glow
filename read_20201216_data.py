import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

infile = 'data/20201216_Lum_Image_Run/20201216_Lum_CellCounts_TOTAL.csv'
cell_line = 'H1048'
drug_name = 'AZD-1152' #'TAK-901'

raw_data = np.genfromtxt(infile, names=True, dtype=None, delimiter=',', encoding='UTF-8')
# print(raw_data)
print(raw_data.dtype.names)

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

t_cell_count = []
cell_count_reps = []
cell_count_mean = []
t_rlu = []
lum_rlu_reps = []
lum_rlu_mean = []
for key in replicates.keys():
    cell_count_reps.append([])
    lum_rlu_reps.append([])
    for rep in replicates[key]:
        cell_count_reps[-1].append(rep['Cell_Count'])
        lum_rlu_reps[-1].append(rep['RLU'])
    t_cell_count.append(replicates[key][0]['TotHour_Image'] / 24.) # time in days
    t_rlu.append(replicates[key][0]['TotHour_Lum'] / 24.)  # time in days
    cell_count_mean.append(np.mean(cell_count_reps[-1], axis=0))
    lum_rlu_mean.append(np.mean(lum_rlu_reps[-1], axis=0))

data_dict = {}
for i in range(len(drug_conc)):
    # time points
    data_dict['t_cell_count_%gnM' % (drug_conc[i]*1e9)] = t_cell_count[i]
    data_dict['t_rlu_%gnM' % (drug_conc[i] * 1e9)] = t_rlu[i]
    # cell counts (mean + replicates)
    data_dict['cells_%gnM_mean' % (drug_conc[i]*1e9)] = cell_count_mean[i]
    for j in range(len(cell_count_reps[i])):
        data_dict['cells_%gnM_rep%d' % (drug_conc[i]*1e9, j)] = cell_count_reps[i][j]
    # luminescence (mean + replicates)
    data_dict['lum_%gnM_mean' % (drug_conc[i]*1e9)] = lum_rlu_mean[i]
    for j in range(len(lum_rlu_reps[i])):
        data_dict['lum_%gnM_rep%d' % (drug_conc[i]*1e9, j)] = lum_rlu_reps[i][j]

print(data_dict.keys())
dose = 0

cmap = cm.get_cmap('jet')
colors = cmap(np.linspace(0,1,len(drug_conc)))

for dose,color in zip(drug_conc*1e9,colors):

    start_idx = 0

    # cell counts
    plt.figure('count')
    init = data_dict['cells_%gnM_mean' % dose][start_idx]
    plt.plot(data_dict['t_cell_count_%gnM' % dose][start_idx:],
             np.log2(data_dict['cells_%gnM_mean' % dose][start_idx:]/init),
             color=color, lw=2, label = '%g nM' % dose)
    # replicates
    # for i in range(3):
    #     init = data_dict['cells_%gnM_rep%d' % (dose,i)][start_idx]
    #     plt.plot(data_dict['t_cell_count_%gnM' % dose], np.log2(data_dict['cells_%gnM_rep%d' % (dose,i)]/init), 'x')
    plt.xlabel('time (day)')
    plt.ylabel('log2(cell count)')
    plt.legend(loc=0, ncol=3)
    plt.tight_layout()

    # luminescence
    plt.figure('lum')
    init = data_dict['lum_%gnM_mean' % dose][start_idx]
    plt.plot(data_dict['t_rlu_%gnM' % dose][start_idx:],
             data_dict['lum_%gnM_mean' % dose][start_idx:]/init,
             color=color, lw=2, label = '%g nM' % dose)
    # replicates
    # for i in range(3):
    #     init = data_dict['lum_%gnM_rep%d' % (dose,i)][start_idx]
    #     plt.plot(data_dict['t_rlu_%gnM' % dose], data_dict['lum_%gnM_rep%d' % (dose,i)]/init, 'x')
    plt.xlabel('time (day)')
    plt.ylabel('luminescence (RLU)')
    plt.legend(loc=0, ncol=3)
    plt.tight_layout()

plt.show()
