import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import Normalize
from scipy.stats import linregress

infile = 'data/20201216_Lum_Image_Run/20201216_Lum_CellCounts_TOTAL_v2.csv'

raw_data = np.genfromtxt(infile, names=True, dtype=None, delimiter=',', encoding='UTF-8')
# print(raw_data)
print(raw_data.dtype.names)

# cell_lines = np.unique([d['Cell_Line'] for d in raw_data])
cell_lines = ['DMS53', 'H1048', 'H841']
# cell_lines = ['DMS53' 'H1048' 'H841']
drug_names = np.unique([d['Drug'] for d in raw_data])
# drug_names = ['barasertib']

t_min = [0]*len(cell_lines)  # days (same length as cell_lines)
t_max = [2, 5, 8]  # [2, 5, 8] # days (same length as cell_lines)

# bad data
# (cell_line, drug_name, drug_conc_in_nM, replicate)
exclude = [('H1048', 'YM-155', 2490, 0),
           ('WM88', True, True, True)]

cmap = cm.get_cmap('coolwarm')

lum_data_out = []
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
        lm = []
        tm = []
        dc = []
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
        #############################################
                # STORE LUMINESCENCE DATA ###
                PROCEED = True
                for x in exclude:
                    if (cell_line == x[0] or x[0] is True) and \
                            (drug_name == x[1] or x[1] is True) and \
                            (drug_conc[i]*1e9 == x[2] or x[2] is True) and \
                            (j == x[3] or x[3] is True):  # first replicate is the bad one
                        PROCEED = False
                if PROCEED:
                    tt = np.array([t for t in t_rlu[i] if t <= t_end])
                    ll = lum_rlu_reps[i][j][:len(tt)]
                    tt = np.array([t for t in tt if t >= t_start])
                    ll = ll[len(ll) - len(tt):]
                    lm.append(ll)
                    tm.append(tt)
                    dc.append(drug_conc[i]*1e9)
        # Store luminescence, time points, and drug concentrations for all cell line/drug combos
        lum_data_out.append((lm, dc, tm, cell_line, drug_name))
        #############################################

        # PLOTS
        vmin = np.floor(np.log10(min(d for d in dc if d > 0)) - 1)
        vmax = np.ceil(np.log10(max(dc)))
        norm = Normalize(vmin=vmin, vmax=vmax)
        fig_cell.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), ax=axs_cell.ravel().tolist(),
                          orientation='vertical', label=r'log$_{10}$(drug conc) (nM)')
        fig_lum.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), ax=axs_lum.ravel().tolist(),
                         orientation='vertical', label=r'log$_{10}$(drug conc) (nM)')

        nconc = 0
        for conc in drug_conc:
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
                # replicates
                for i in range(n_replicates):
                    #####
                    t_count = data_dict['t_cell_count_%gnM' % conc]
                    count = data_dict['cells_%gnM_rep%d' % (conc, i)]
                    # t_count = np.array([t for t in data_dict['t_cell_count_%gnM' % conc] if t <= t_end])
                    # count = data_dict['cells_%gnM_rep%d' % (conc, i)][:len(t_count)]
                    # t_count = np.array([t for t in t_count if t >= t_start])
                    # count = count[len(count) - len(t_count):]
                    #####
                    label = '%g nM' % conc if i == 0 else None
                    with np.errstate(divide='ignore'):  # suppress the divide by zero warning
                        log10_conc = max(np.log10(conc), vmin)
                        color = cmap((log10_conc - vmin) / (vmax - vmin))
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
                    ax_cell.set_ylabel('doublings')
                    # ax_cell.set_ylabel('cell count')
                    ax_dip.set_ylabel('DIP rate (1/day)')
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
                    with np.errstate(divide='ignore'):
                        log10_conc = max(np.log10(conc), vmin)
                        color = cmap((log10_conc - vmin) / (vmax - vmin))
                    ax_lum.plot(t_lum, lum, color=color, lw=1, label=label)
                    # ax_lum.plot(data_dict['t_rlu_%gnM' % conc],
                    #             data_dict['lum_%gnM_rep%d' % (conc, i)],
                    #             color=color, lw=1, label=label)
                if row == nrows - 1 or (row + 1) * (col + 1) + ncols > len(drug_names):
                    ax_lum.set_xlabel('time (day)')
                if col == 0:
                    ax_lum.set_ylabel('luminescence (RLU)')
                # ax_lum.legend(loc=0, ncol=3)
                ax_lum.set_xlim(xmin=t_start, xmax=t_end)

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

##############
# lum_data_out is the object we want to send to PSO, but the data needs to be cleaned up first
# Elements of lum_data_out are (lum, drug_conc, time, cell_line, drug_name)
##############
pso_rt_glo = __import__('PSO-RT-Glo')

cmap = cm.get_cmap('coolwarm')

cell_lines = np.unique([data[3] for data in lum_data_out])
drug_names = np.unique([data[4] for data in lum_data_out])

for cell_line in cell_lines:
    ncols = int(np.sqrt(len(drug_names)))
    nrows = int(np.ceil(len(drug_names) / ncols))
    fig_lum_sim, axs_lum_sim = plt.subplots(nrows=nrows, ncols=ncols, sharex='col', sharey='row', figsize=(7, 7),
                                            constrained_layout=True, num='%s luminescence (model fits)' % cell_line)
    fig_cell_sim, axs_cell_sim = plt.subplots(nrows=nrows, ncols=ncols, sharex='col', sharey='row', figsize=(7, 7),
                                              constrained_layout=True, num='%s cell counts (model predicted)' %
                                                                           cell_line)
    fig_lum_exp, axs_lum_exp = plt.subplots(nrows=nrows, ncols=ncols, sharex='col', sharey='row', figsize=(7, 7),
                                            constrained_layout=True, num='%s luminescence (experiment)' % cell_line)
    if ncols == 1 and nrows == 1:
        axs_lum_sim = np.array(axs_lum_sim)
        axs_lum_exp = np.array(axs_lum_exp)
    if ncols == 1:
        axs_lum_sim.shape = (nrows, 1)
        axs_lum_exp.shape = (nrows, 1)
    for n, drug_name in enumerate(drug_names):
        print(cell_line, drug_name)
        exp_data = [d for d in lum_data_out if d[3] == cell_line and d[4] == drug_name][0]
        lum_exp = exp_data[0]
        drug_conc = exp_data[1]
        time = exp_data[2]
        dc, idx = np.unique(drug_conc, return_index=True)
        lum_exp_avg = np.array([np.mean(lum_exp[idx[i]:idx[i + 1]], axis=0) for i in range(len(idx) - 1)] +
                               [np.mean(lum_exp[idx[-1]:], axis=0)])
        lum_exp_var = np.array([np.var(lum_exp[idx[i]:idx[i + 1]], axis=0) for i in range(len(idx) - 1)] +
                               [np.var(lum_exp[idx[-1]:], axis=0)])
        tm = [time[i] for i in idx]
        # run PSO
        m, nC0, a, b, c = pso_rt_glo.run_pso(lum_exp_avg, dc, tm, 1 / lum_exp_var)
        print(m, nC0, a, b, c)  # m, nC0, kdiv-kdeath, kdiv*-kdeath*, koff/kon
        # plot simulated and experimental data
        col = n % ncols
        row = int(n / ncols)
        vmin = np.floor(np.log10(min(d for d in dc if d > 0)) - 1)
        vmax = np.log10(max(dc))
        norm = Normalize(vmin=vmin, vmax=vmax)
        # simulated luminescence
        ax_lum_sim = axs_lum_sim[row, col]
        fig_lum_sim.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), ax=axs_lum_sim.ravel().tolist(),
                             orientation='vertical', label=r'log$_{10}$(drug conc) (nM)')
        ax_lum_sim.set_title('%s+%s' % (cell_line, drug_name))
        ax_lum_sim.set_xlabel('time (day)')
        ax_lum_sim.set_ylabel('luminescence (RLU)')
        # simulated cell counts
        ax_cell_sim = axs_cell_sim[row, col]
        fig_cell_sim.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), ax=axs_cell_sim.ravel().tolist(),
                             orientation='vertical', label=r'log$_{10}$(drug conc) (nM)')
        ax_cell_sim.set_title('%s+%s' % (cell_line, drug_name))
        ax_cell_sim.set_xlabel('time (day)')
        ax_cell_sim.set_ylabel('doublings')
        # experimental luminescence
        ax_lum_exp = axs_lum_exp[row, col]
        fig_lum_exp.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), ax=axs_lum_exp.ravel().tolist(),
                             orientation='vertical', label=r'log$_{10}$(drug conc) (nM)')
        ax_lum_exp.set_title('%s+%s' % (cell_line, drug_name))
        ax_lum_exp.set_xlabel('time (day)')
        ax_lum_exp.set_ylabel('luminescence (RLU)')
        #####
        for i in range(len(dc)):
            # plot simulated luminescence
            lum_sim = np.maximum(
                nC0 * (1 - np.exp(a * tm[i])) + m * tm[i] + nC0 *
                np.exp(((c * a + dc[i] * b) / (c + dc[i])) * tm[i]),
                np.zeros(len(tm[i]))
            )
            log10_conc = max(np.log10(dc[i]), vmin)
            color = cmap((log10_conc - vmin) / (vmax - vmin))
            ax_lum_sim.plot(tm[i], lum_sim, color=color, lw=1)
            #####
            # plot simulated cell population doublings
            t_cell = np.linspace(0, 8, 101)
            ln_cell_ratio = ((c * a + dc[i] * b) / (c + dc[i])) * t_cell  #tm[i]
            log2_cell_ratio = ln_cell_ratio / np.log(2)
            # ax_cell_sim.plot(tm[i], log2_cell_ratio, color=color, lw=1)
            ax_cell_sim.plot(t_cell, log2_cell_ratio, color=color, lw=1)
            #####
            # plot experimental luminescence
            idx_end = len(lum_exp) if i == len(dc)-1 else idx[i+1]
            for j in range(idx[i], idx_end):
                ax_lum_exp.plot(time[j], lum_exp[j], color=color, lw=1)
        # adjust y-axis limits
        ylim = ax_lum_exp.get_ylim()
        ymin = 0 if ylim[0] > 0 else ylim[0]
        ax_lum_sim.set_ylim(ymin, ylim[1])
        ax_lum_exp.set_ylim(bottom=ymin)

plt.show()


quit()

#################
# OLD CODE BELOW
#################

for data in lum_data_out[:]:  # (lm, dc, tm, cell_line, drug_name)
    lum, drug_conc, time, cell_line, drug_name = data
    dc, idx = np.unique(drug_conc, return_index=True)
    lum_avg = np.array([np.mean(lum[idx[i]:idx[i + 1]], axis=0) for i in range(len(idx) - 1)] +
                       [np.mean(lum[idx[-1]:], axis=0)])
    lum_var = np.array([np.var(lum[idx[i]:idx[i + 1]], axis=0) for i in range(len(idx) - 1)] +
                       [np.var(lum[idx[-1]:], axis=0)])
    tm = [time[i] for i in idx]

    print('%s+%s' % (cell_line, drug_name))
    colors = cmap(np.linspace(0, 1, len(dc)))
    norm = Normalize(vmin=dc[0] / 1e3, vmax=dc[-1] / 1e3)
    #####
    # plot experimental data
    plt.figure('%s+%s: Experiment' % (cell_line, drug_name))
    plt.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), ax=plt.axes(),
                 orientation='vertical', label=r'drug conc ($\mu$M)')
    plt.xlabel('time (day)')
    plt.ylabel('luminescence (RLU)')
    #####
    # run PSO and plot simulated data
    m, nC0, a, b, c = pso_rt_glo.run_pso(lum_avg, dc, tm, 1/lum_var)
    print(m, nC0, a, b, c)  # m, nC0, kdiv-kdeath, kdiv*-kdeath*, koff/kon
    plt.figure('%s+%s: Simulation' % (cell_line, drug_name))
    plt.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), ax=plt.axes(),
                 orientation='vertical', label=r'drug conc ($\mu$M)')
    plt.xlabel('time (day)')
    plt.ylabel('luminescence (RLU)')
    #####
    xmin = min(tm[0])
    xmax = max(tm[0])
    # xmin = min(time[0])
    # xmax = max(time[0])
    for i in range(len(dc)):
        plt.figure('%s+%s: Experiment' % (cell_line, drug_name))
        # plot mean luminescence
        plt.plot(tm[i], lum_avg[i], color=colors[i])
        # plot replicates
        if i < len(dc)-1:
            plt.plot(time[idx[i]:idx[i+1]], lum[idx[i]:idx[i+1]], 'o', color=colors[i])
        else:
            plt.plot(time[idx[i]:], lum[idx[i]:], 'o', color=colors[i])
        if min(tm[i]) < xmin:
            xmin = min(tm[i])
        if max(tm[i]) > xmax:
            xmax = max(tm[i])
        ##########
        # plot weights (inverse of variances)
        plt.figure('%s+%s: Weights' % (cell_line, drug_name))
        plt.plot(tm[i], np.log10(1 / lum_var[i]), '-o', lw=2, color=colors[i])
        plt.xlabel('time')
        plt.ylabel('weight (1/var)')
        plt.tight_layout()
        ##########
        plt.figure('%s+%s: Simulation' % (cell_line, drug_name))
        lum_sim = np.maximum(
            nC0 * (1 - np.exp(a * tm[i])) + m * tm[i] + nC0 *
            np.exp(((c * a + dc[i] * b) / (c + dc[i])) * tm[i]),
            np.zeros(len(tm[i]))
        )
        plt.plot(tm[i], lum_sim, color=colors[i])
    # set x-axis limits
    plt.figure('%s+%s: Experiment' % (cell_line, drug_name))
    xlim = plt.xlim((np.floor(xmin), np.ceil(xmax)))
    ylim = plt.ylim()
    plt.figure('%s+%s: Simulation' % (cell_line, drug_name))
    plt.xlim(xlim)
    plt.ylim(ylim)
#############

plt.show()
