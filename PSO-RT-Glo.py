from pysb import *
from pysb.integrate import Solver
from pysb.simulator import ScipyOdeSimulator
import matplotlib.pyplot as plt
import numpy as np
import itertools
import sympy
from sklearn.linear_model import LinearRegression
from scipy.stats import linregress
import math
from simplepso.pso import PSO

# The Below code applies only to one cell line
a = 0.5  # Kdiv - kdeath
b = -0.584  # Kdiv* - kdeath* orig:-0.584
c = 2.117  # Koff/kon
m = 2500.0 / 7
nC0 = 500
# #print(t)
np.random.seed(13)
# Getting the m and nC0 Parameter Values
noise1 = np.random.normal(0, 125, 25)
noise2 = np.random.normal(0, 125, 25)
noise3 = np.random.normal(0, 125, 25)
t1 = np.linspace(0, 7, 25)
t2 = np.linspace(0, 7, 25)
t3 = np.linspace(0, 7, 25)
conc1 = 0
conc2 = 1
conc3 = 10
c_count1 = 50 * np.exp(((c * a + conc1 * b) / (c + conc1)) * t1)
c_count2 = 50 * np.exp(((c * a + conc2 * b) / (c + conc2)) * t2)
c_count3 = 50 * np.exp(((c * a + conc3 * b) / (c + conc3)) * t3)
lum1 = np.maximum(nC0 * (1 - np.exp(a * t1)) + m * t1 + nC0 / 50 * c_count1 + noise1, np.zeros(len(t1)))
lum2 = np.maximum(nC0 * (1 - np.exp(a * t2)) + m * t2 + nC0 / 50 * c_count2 + noise2, np.zeros(len(t2)))
lum3 = np.maximum(nC0 * (1 - np.exp(a * t3)) + m * t3 + nC0 / 50 * c_count3 + noise3, np.zeros(len(t3)))

lum1_prime = nC0 / 50 * c_count1
lum2_prime = nC0 / 50 * c_count2
lum3_prime = nC0 / 50 * c_count3
plt.plot(t1, lum1)
# plt.plot(t1, lum1_prime)

plt.plot(t2, lum2)
# plt.plot(t2, lum2_prime)
plt.plot(t3, lum3)
# plt.plot(t2, lum3_prime)

plt.figure()
plt.plot(t1, np.log2(c_count1), label='Conc = %d' % conc1)
plt.plot(t2, np.log2(c_count2), label='Conc = %d' % conc2)
plt.plot(t3, np.log2(c_count3), label='Conc = %d' % conc3)
plt.legend(loc=0)
plt.show()

# lum1 = 500 + (2500/7)*t1 + noise1
# lum2 = 500 + (1500/7)*t2 + noise2
# lum3 = 500 + (500/7)*t3 + noise3

###################################
time = np.array([t1, t2, t3])
drug_conc = np.array([conc1, conc2, conc3])
lum = np.array([lum1, lum2, lum3])
###################################
# #print(data)
# fit1 = linregress(t, lum1)
# fit2 = linregress(t, lum2)
# fit3 = linregress(t, lum3)

# y_pred1 = fit1.intercept+fit1.slope*t
# y_pred2 = fit2.intercept+fit2.slope*t
# y_pred3 = fit3.intercept+fit3.slope*t
# plt.plot(t, lum1, "+r", label = 'Synthetic Data')
# plt.plot(t, lum2, "+r")
# plt.plot(t, lum3, "+r")
# plt.plot(t, y_pred1, "b", label = 'Linear Fit')
# plt.plot(t, y_pred2, "b")
# plt.plot(t, y_pred3, "b")
# plt.xlabel('Time (Days)')
# plt.ylabel('Luminescence')
# plt.legend(loc = 0)
# plt.show()
# #print(noise)

# # Experimental Data


# The following PSO Code only takes in 1 slice (Drug and Cell Combo)
# For the Fitting Equation, a = kdiv-death, b = kdiv*-kdeath*, and c = koff/on
# Therefore the fitting equation would look like: l = nC0*(1-e**(a*t)) + mt + nC0*2.71**(((c*a + d*b)/(c+d))*t)

SSD_list = []
a_list = []
b_list = []
c_list = []


def cost(params):
    a = params[0]
    b = params[1]
    c = params[2]

    # Using Model Equation
    SSD = np.empty(len(drug_conc), dtype=object)

    for i, d in enumerate(drug_conc):
        traj = np.maximum(
            nC0 * (1 - np.exp(a * time[i])) + m * time[i] + nC0 * np.exp(((c * a + d * b) / (c + d)) * time[i]),
            np.zeros(len(time[i])))
        SSD[i] = sum((lum[i] - traj) ** 2)

    SSD_list.append(sum(SSD))
    a_list.append(a)
    b_list.append(b)
    c_list.append(c)

    return sum(SSD)


# Make array of drug conc
def run_pso(lum, drug_conc, time):
    lum = np.array(lum)
    drug_conc = np.array(drug_conc)
    time = np.array(time)

    global m
    m = 0
    global nC0
    nC0 = 0
    indices = np.where(drug_conc == 0)[0]
    for idx in indices:
        fit = linregress(time[idx], lum[idx])
        m += fit.slope
        nC0 += fit.intercept
    m /= len(indices)
    nC0 /= len(indices)

    print(m)
    print(nC0)

    pso = PSO(save_sampled=True, verbose=True, shrink_steps=False)
    pso.set_start_position([2, 2, 2])

    # allows particles to move +/- 2 orders of magnitude
    pso.set_bounds(2)
    # sets maximum speed that a particle can travel
    pso.set_speed(-.1, .1)
    pso.run(
        num_particles=1,
        num_iterations=1000,
        stop_threshold=0,
        num_processors=1,
        max_iter_no_improv=100,
        cost_function=cost
    )


##################################################
# Main Code

##### Input data #####

# lum = [[6,10,13,14],
#         [7,11,14],
#         [1,3,5,7,9],
#         [1,2,4,9,10],
#         [0.5,1]]

# drug_conc = [0,0,0.1,0.1,0.25]

# time = [[1,2,3,4],
#         [1,2,3],
#         [1,2,3,4,5],
#         [1,2,3,4],
#         [1,2]]

run_pso(lum, drug_conc, time)
#####################


plt.plot(range(len(SSD_list)), np.log10(np.array(SSD_list)))
plt.xlabel('Iteration')
plt.ylabel('log10(Cost)')
plt.figure()
plt.plot(range(len(a_list)), a_list, ".", label="Kdiv-Kdeath")
plt.plot(range(len(b_list)), b_list, ".", label="Kdiv*-Kdeath*")
plt.plot(range(len(c_list)), c_list, ".", label="Koff/Kon")
plt.xlabel('Iteration')
plt.ylabel('Param Value')
plt.legend(loc=0)
print(np.log10(SSD_list[-1]))
print('Final Values:')
for i in range(-5, 0):
    print(a_list[i], b_list[i], c_list[i])
print("Initial Values")
print(a, b, c)

plt.figure()
colors = ['r', 'b', 'g']
a = a_list[-1]
b = b_list[-1]
c = c_list[-1]
for i, d in enumerate(drug_conc):
    plt.plot(time[i], lum[i], '+', color=colors[i])
    lum_model = np.maximum(
        nC0 * (1 - np.exp(a * time[i])) + m * time[i] + nC0 * np.exp(((c * a + d * b) / (c + d)) * time[i]),
        np.zeros(len(time[i])))
    plt.plot(time[i], lum_model, color=colors[i])
    plt.xlabel('Time')
    plt.ylabel('Luminescence')

plt.show()
