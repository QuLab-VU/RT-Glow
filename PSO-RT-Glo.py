import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import linregress
from simplepso.pso import PSO

sse_list = []
a_list = []
b_list = []
c_list = []


def cost(params):

    a = params[0]
    b = params[1]
    c = params[2]

    sse = np.empty(len(_drug_conc), dtype=object)

    for i, d in enumerate(_drug_conc):
        traj = np.maximum(
            _nC0 * (1 - np.exp(a * _time[i])) + _m * _time[i] + _nC0 * np.exp(((c * a + d * b) / (c + d)) * _time[i]),
            np.zeros(len(_time[i])))
        sse[i] = sum((_lum[i] - traj) ** 2)

    #####
    sse_list.append(sum(sse))
    a_list.append(a)
    b_list.append(b)
    c_list.append(c)
    #####

    return sum(sse)


def run_pso(lum, drug_conc, time):

    # create global versions of lum, drug_conc, and time variables that can be used in cost function
    global _lum, _drug_conc, _time
    _lum = np.array(lum)
    _drug_conc = np.array(drug_conc)
    _time = np.array(time)

    # estimate slope m and y-intercept nC0 from control luminescence data
    global _m, _nC0
    _m = 0
    _nC0 = 0
    # find indices of control experiments
    indices = np.where(_drug_conc == 0)[0]
    for idx in indices:
        fit = linregress(_time[idx], _lum[idx])
        _m += fit.slope
        _nC0 += fit.intercept
    _m /= len(indices)
    _nC0 /= len(indices)

    # run PSO to estimate the other three parameters: kdiv-kdeath, kdiv*-kdeath*, koff/kon
    pso = PSO(save_sampled=True, verbose=True, shrink_steps=False)
    pso.set_start_position([2, 2, 2])
    # allows particles to move +/- 2 orders of magnitude
    pso.set_bounds(2)
    # sets maximum speed that a particle can travel
    pso.set_speed(-0.1, 0.1)
    pso.run(
        num_particles=1,
        num_iterations=1000,
        stop_threshold=0,
        num_processors=1,
        max_iter_no_improv=100,
        cost_function=cost
    )


##################################################

def run_example():

    # EXAMPLE: CREATE SOME SYNTHETIC DATA

    # model parameters
    a = 0.7  # 0.5  # Kdiv - kdeath
    b = 0.1  # -0.584  # Kdiv* - kdeath* orig:-0.584
    c = 2.9  # 2.117  # Koff/kon
    m = 2500.0 / 7  # slope of control luminescence time course
    nC0 = 500  # initial luminescence

    # drug concentrations
    conc1 = 0
    conc2 = 1
    conc3 = 10

    # time points for different drug concentations
    t1 = np.linspace(0, 7, 25)
    t2 = np.linspace(0, 7, 25)
    t3 = np.linspace(0, 7, 25)

    # cell counts for different drug concentrations
    c_count1 = 50 * np.exp(((c * a + conc1 * b) / (c + conc1)) * t1)
    c_count2 = 50 * np.exp(((c * a + conc2 * b) / (c + conc2)) * t2)
    c_count3 = 50 * np.exp(((c * a + conc3 * b) / (c + conc3)) * t3)

    # add some noise
    np.random.seed(13)
    noise1 = np.random.normal(0, 125, 25)
    noise2 = np.random.normal(0, 125, 25)
    noise3 = np.random.normal(0, 125, 25)

    # luminescence trajectories for different drug concentrations
    lum1 = np.maximum(nC0 * (1 - np.exp(a * t1)) + m * t1 + nC0 / 50 * c_count1 + noise1, np.zeros(len(t1)))
    lum2 = np.maximum(nC0 * (1 - np.exp(a * t2)) + m * t2 + nC0 / 50 * c_count2 + noise2, np.zeros(len(t2)))
    lum3 = np.maximum(nC0 * (1 - np.exp(a * t3)) + m * t3 + nC0 / 50 * c_count3 + noise3, np.zeros(len(t3)))

    # plot synthetic data
    fig, axs = plt.subplots(1, 2, figsize=(12.8, 4.8))
    fig.suptitle('Synthetic Data')

    axs[0].plot(t1, lum1, label='Conc = %d' % conc1)
    axs[0].plot(t2, lum2, label='Conc = %d' % conc2)
    axs[0].plot(t3, lum3, label='Conc = %d' % conc3)
    axs[0].set_xlabel('time')
    axs[0].set_ylabel('luminescence')
    axs[0].legend(loc=0)

    axs[1].plot(t1, np.log2(c_count1), label='Conc = %d' % conc1)
    axs[1].plot(t2, np.log2(c_count2), label='Conc = %d' % conc2)
    axs[1].plot(t3, np.log2(c_count3), label='Conc = %d' % conc3)
    axs[1].set_xlabel('time')
    axs[1].set_ylabel('log2(cell count)')
    axs[1].legend(loc=0)

    # data objects to be input to PSO code
    time = np.array([t1, t2, t3])
    drug_conc = np.array([conc1, conc2, conc3])
    lum = np.array([lum1, lum2, lum3])

    # run PSO
    run_pso(lum, drug_conc, time)

    # cost vs. iteration
    plt.figure()
    plt.plot(range(len(sse_list)), np.log10(np.array(sse_list)))
    plt.xlabel('Iteration')
    plt.ylabel('log10(Cost)')

    # parameter values vs. iteration
    plt.figure()
    plt.plot(range(len(a_list)), a_list, ".", label="Kdiv-Kdeath")
    plt.plot(range(len(b_list)), b_list, ".", label="Kdiv*-Kdeath*")
    plt.plot(range(len(c_list)), c_list, ".", label="Koff/Kon")
    plt.xlabel('Iteration')
    plt.ylabel('Param Value')
    plt.legend(loc=0)

    print(np.log10(sse_list[-1]))
    print('Final Values:')
    for i in range(-5, 0):
        print(a_list[i], b_list[i], c_list[i])
    print("Initial Values")
    print(a, b, c)

    # luminescence fits and inferred cell counts
    fig, axs = plt.subplots(1, 2, figsize=(12.8, 4.8))
    fig.suptitle('Luminescence Fits and Inferred Cell Counts')
    colors = ['r', 'b', 'g']
    a = a_list[-1]
    b = b_list[-1]
    c = c_list[-1]
    for i, d in enumerate(drug_conc):
        axs[0].plot(time[i], lum[i], '+', color=colors[i])
        cells_model = 50 * np.exp(((c * a + d * b) / (c + d)) * time[i])
        lum_model = np.maximum(
            nC0 * (1 - np.exp(a * time[i])) + m * time[i] + nC0/50 * cells_model, np.zeros(len(time[i])))
        axs[0].plot(time[i], lum_model, color=colors[i])
        axs[1].plot(time[i], np.log2(cells_model), color=colors[i])
        plt.xlabel('Time')
        plt.ylabel('Luminescence')

    plt.show()


run_example()
