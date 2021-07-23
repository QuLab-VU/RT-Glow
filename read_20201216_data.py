import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import Normalize
from scipy.stats import linregress

# x = np.array([20, 10])
# y = np.array([4,2])
# print(x/y)
# print(x+y)
# print(x*y)
# quit()

infile = 'data/20201216_Lum_Image_Run/20201216_Lum_CellCounts_TOTAL_v2.csv'

raw_data = np.genfromtxt(infile, names=True, dtype=None, delimiter=',', encoding='UTF-8')
# print(raw_data)
print(raw_data.dtype.names)

# cell_lines = np.unique([d['Cell_Line'] for d in raw_data])
cell_lines = ['DMS53', 'H1048', 'H841']
# cell_lines = ['DMS53' 'H1048' 'H841']
drug_names = np.unique([d['Drug'] for d in raw_data])
# drug_names = ['barasertib']

t_min = [1]*len(cell_lines) # days (same length as cell_lines)
t_max = [8, 8, 8]  # [2, 5, 8] # days (same length as cell_lines)

# bad data
exclude = [('H1048', 'YM-155', 2490),
           ('WM88', True, True)]

cmap = cm.get_cmap('coolwarm')

for cell_line, t_start, t_end in zip(cell_lines, t_min, t_max):

    ncols = int(np.sqrt(len(drug_names)))
    nrows = int(np.ceil(len(drug_names) / ncols))

    fig_cell, axs_cell = plt.subplots(nrows=nrows, ncols=ncols, sharex='col', sharey='row', figsize=(7, 7),
                                      constrained_layout=True, num='%s cell counts' % cell_line)

    fig_dip, axs_dip = plt.subplots(nrows=nrows, ncols=ncols, sharex='col', sharey='row', figsize=(7, 7),
                                      constrained_layout=True, num='%s DIP rates' % cell_line)

    fig_lum, axs_lum = plt.subplots(nrows=nrows, ncols=ncols, sharex='col', sharey='row', figsize=(7, 7),
                                    constrained_layout=True, num='%s luminescence' % cell_line)

    if ncols == 1 and nrows == 1:
        axs_cell = np.array(axs_cell)
        axs_dip = np.array(axs_dip)
        axs_lum = np.array(axs_lum)
    if ncols == 1:
        axs_cell.shape = (nrows, 1)
        axs_dip.shape = (nrows, 1)
        axs_lum.shape = (nrows, 1)

    # delete unused plots
    for ax in axs_cell.flat[len(drug_names):]:
        fig_cell.delaxes(ax)
    for ax in axs_dip.flat[len(drug_names):]:
        fig_dip.delaxes(ax)
    for ax in axs_lum.flat[len(drug_names):]:
        fig_lum.delaxes(ax)

    ndrugs = 0
    for n, drug_name in enumerate(drug_names):
        print(cell_line, drug_name)

        row = int(n/ncols)
        col = n % ncols
        ax_cell = axs_cell[row, col]
        ax_dip = axs_dip[row, col]
        ax_lum = axs_lum[row, col]

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
            t_cell_count.append(replicates[key][0]['TotHour_Image'] / 24.)  # time in days
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

        # PLOTS
        colors = cmap(np.linspace(0, 1, len(drug_conc)))
        norm = Normalize(vmin=drug_conc[0]*1e6, vmax=drug_conc[-1]*1e6)
        fig_cell.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), ax=axs_cell.ravel().tolist(),
                          orientation='vertical', label=r'drug conc ($\mu$M)')
        fig_lum.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), ax=axs_lum.ravel().tolist(),
                         orientation='vertical', label=r'drug conc ($\mu$M)')

        nconc = 0
        for conc, color in zip(drug_conc, colors):
            n_replicates = len(replicates[conc])
            conc = conc*1e9
            # check for bad data
            PROCEED = True
            for x in exclude:
                if (cell_line == x[0] or x[0] is True) and \
                        (drug_name == x[1] or x[1] is True) and \
                        (conc == x[2] or x[2] is True):
                    PROCEED = False
            if PROCEED:
                nconc += 1
                # cell counts
                ax_cell.set_title('%s+%s' % (cell_line, drug_name))
                ax_dip.set_title('%s+%s' % (cell_line, drug_name))
                # mean
                # init = data_dict['cells_%gnM_mean' % conc][start_idx]
                # plt.plot(data_dict['t_cell_count_%gnM' % conc][start_idx:],
                #          np.log2(data_dict['cells_%gnM_mean' % conc][start_idx:]/init),
                #          color=color, lw=1, label = '%g nM' % conc)
                # replicates
                for i in range(n_replicates):
                    #####
                    t_count = np.array([t for t in data_dict['t_cell_count_%gnM' % conc] if t <= t_end])
                    count = data_dict['cells_%gnM_rep%d' % (conc, i)][:len(t_count)]
                    t_count = np.array([t for t in t_count if t >= t_start])
                    count = count[len(count) - len(t_count):]
                    #####
                    label = '%g nM' % conc if i == 0 else None
                    ax_cell.plot(t_count, np.log2(count / count[0]), color=color, lw=1, label=label)
                    # ax_cell.plot(data_dict['t_cell_count_%gnM' % conc],
                    #              np.log2(data_dict['cells_%gnM_rep%d' % (conc, i)]/init),
                    #              color=color, lw=1, label=label)
                    # ax_cell.plot(data_dict['t_cell_count_%gnM' % conc],
                    #              data_dict['cells_%gnM_rep%d' % (conc, i)],
                    #              color=color, lw=1, label=label)
                    #####
                    # calculate DIP rate
                    fit = linregress(t_count, np.log2(count / count[0]))
                    # ax_dip.plot(t_count, fit.slope * (t_count - t_count[0]) + 1, color=color, lw=1)
                    # print(np.log10(conc/1e9))
                    if conc == 0.0:
                        d = np.min([c for c in drug_conc if c > 0.0]) / 100
                    else:
                        d = conc/1e9
                    ax_dip.plot(np.log10(d), fit.slope, 'kx', mfc='None')
                    ax_dip.set_ylim(ymin=-0.5, ymax=1)
                if row == nrows - 1 or (row + 1) * (col + 1) + ncols > len(drug_names):
                    ax_cell.set_xlabel('time (day)')
                    ax_dip.set_xlabel('log10[conc (M)]')
                if col == 0:
                    ax_cell.set_ylabel('log2(cell count)')
                    # ax_cell.set_ylabel('cell count')
                    ax_dip.set_ylabel('DIP rate')
                # ax_cell.legend(loc=0, ncol=3)

                # luminescence
                # plt.figure('lum %s+%s' % (cell_line, drug_name))
                ax_lum.set_title('%s+%s' % (cell_line, drug_name))
                # mean
                # init = data_dict['lum_%gnM_mean' % conc][start_idx]
                # plt.plot(data_dict['t_rlu_%gnM' % conc][start_idx:],
                #          data_dict['lum_%gnM_mean' % conc][start_idx:]/init,
                #          color=color, lw=1, label = '%g nM' % conc)
                # replicates
                for i in range(n_replicates):
                    label = '%g nM' % conc if i == 0 else None
                    #####
                    t_lum = np.array([t for t in data_dict['t_rlu_%gnM' % conc] if t <= t_end])
                    lum = data_dict['lum_%gnM_rep%d' % (conc, i)][:len(t_lum)]
                    t_lum = np.array([t for t in t_lum if t >= t_start])
                    lum = lum[len(lum) - len(t_lum):]
                    #####
                    ax_lum.plot(t_lum, lum, color=color, lw=1, label=label)
                    # ax_lum.plot(data_dict['t_rlu_%gnM' % conc],
                    #             data_dict['lum_%gnM_rep%d' % (conc, i)],
                    #             color=color, lw=1, label=label)
                if row == nrows - 1 or (row + 1) * (col + 1) + ncols > len(drug_names):
                    ax_lum.set_xlabel('time (day)')
                if col == 0:
                    ax_lum.set_ylabel('luminescence (RLU)')
                # ax_lum.legend(loc=0, ncol=3)

        # if there's no good concentration data, delete the plot
        if nconc == 0:
            fig_cell.delaxes(ax_cell)
            fig_lum.delaxes(ax_lum)
        else:
            ndrugs += 1

    # if there are no plots, delete the figure
    if ndrugs == 0:
        fig_cell.clf()
        fig_lum.clf()
        plt.close(fig_cell)
        plt.close(fig_lum)

plt.show()
