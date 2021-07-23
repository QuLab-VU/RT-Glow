import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import Normalize

# Cu -> Cu + Cu, kdiv_u
# Cu -> 0, kdth_u
# Cd -> Cd + Cd, kdiv_d
# Cd -> 0, kdth_d
# Cu + D -> Cd + D, kon
# Cd -> Cu, koff

# ln(CT/CT_0) = ( ( (koff/kon)*(kdiv_u-kdth_u) + [D]*(kdiv_d-kdth_d) ) * t ) / ( koff/kon + [D] )

kdiv_u = 0.065*np.log(2)
kdth_u = 0.005*np.log(2)
kdiv_d = 0.01*np.log(2)
kdth_d = 0.04*np.log(2)
kon = 1000
koff = 10
CT_0 = 1

conc = 10**np.arange(-6, 0.1, 0.1)
t = np.linspace(0, 120, 1201)

n = 50
m = 100

cmap = cm.get_cmap('coolwarm')
colors = cmap(np.linspace(0, 1, len(conc)))
norm = Normalize(vmin=conc[0], vmax=conc[-1])

fig, axs = plt.subplots(nrows=1, ncols=2, constrained_layout=True, figsize=(10, 5))
cb = fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), ax=axs.ravel().tolist(),
                  orientation='vertical', label=r'drug conc ($\mu$M)')
cb.ax.tick_params(labelsize=14)
cb.ax.yaxis.label.set_size(16)

for color, D in zip(colors, conc):
    print(D)
    CT = CT_0 * np.exp((((koff/kon)*(kdiv_u-kdth_u) + D*(kdiv_d-kdth_d)) * t) / (koff/kon + D))
    # plt.figure('cells')
    axs[0].plot(t/24., np.log2(CT/CT_0), lw=2, color=color)
    # plt.plot(t, np.log2(CT), lw=2)
    axs[0].set_xlabel('time (day)', fontsize=16)
    axs[0].set_ylabel('cell count', fontsize=16)
    axs[0].tick_params(axis='both', which='major', labelsize=16)
    #####
    lum = n*CT - (n*CT_0*(np.exp((kdiv_u-kdth_u)*t) - 1) - m*t)
    # lum = n*CT_0 + m*t
    # plt.figure('lum')
    axs[1].plot(t/24., lum/lum[0], lw=2, color=color)
    axs[1].set_xlabel('time(day)', fontsize=16)
    axs[1].set_ylabel('luminescence', fontsize=16)
    axs[1].tick_params(axis='both', which='major', labelsize=16)

plt.show()


