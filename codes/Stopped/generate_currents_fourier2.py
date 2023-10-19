import matplotlib.pyplot as plt
from math import pi
import numpy as np
from scipy.optimize import curve_fit

import sys
sys.path.append('../') 
from utils import get_data

path = "stopped_data/"
filename = "stopped_data_test2.txt"
times,positions, velocities, currents, uqs, supplys, errors_pos = get_data(filename, path)

# Consolidate Current of bigger magnitude
def get_consolidate_data(positions, currents):
    positions_temp, currents_temp = [], []
    for i in range(len(positions)):
        if positions[i] > 0:
            pos = abs(positions[i])%(2*pi)
        else:
            pos = (2*pi)-(abs(positions[i])%(2*pi))

        if pos not in positions_temp:
            positions_temp.append(pos)
            currents_temp.append(currents[i])
        else:
            ind = positions_temp.index(pos)
            if abs(currents[i]) > abs(currents_temp[ind]):
                currents_temp[ind] = currents[i]

    return positions_temp, currents_temp

index_max = positions.index(max(positions))
positions_reverse = positions[index_max:]
currents_reverse = currents[index_max:]
positions = positions[ :index_max + 1]
currents = currents[ :index_max + 1]

positions, currents = get_consolidate_data(positions, currents)
positions_reverse, currents_reverse = get_consolidate_data(positions_reverse, currents_reverse)

plt.plot(positions, currents, ".", label = "Forward Data")
plt.plot(positions_reverse, currents_reverse, ".", label = "Reverse Data")
plt.title("Cogging mapping with motor stopped")
plt.xlabel("Position [rad]")
plt.ylabel("Current [A]")
plt.legend()
plt.show()

positions += positions_reverse
currents += currents_reverse

# Combine os arrays em pares (posição, corrente)
combined = list(zip(positions, currents))
sorted_combined = sorted(combined, key=lambda x: x[0])
positions, currents = zip(*sorted_combined)
positions = list(positions)
currents = list(currents)

# Get currents uniform in position
INTERVAL_POSITION = 0.002

NUM_POINTS = int(2*pi / INTERVAL_POSITION)
uniform_positions_desired = np.linspace(0, 2*pi, NUM_POINTS)
uniform_currents = []
uniform_positions = []

for i in range(len(uniform_positions_desired)):
    temp_cs = []
    left_lim = uniform_positions_desired[i-1]  
    if i != len(uniform_positions_desired) - 1:
        right_lim = uniform_positions_desired[i+1]        
    else:
        right_lim = uniform_positions_desired[1]

    if i == 0 or i == len(uniform_positions_desired) - 1:
        for j in range(len(positions)):
            if positions[j] > left_lim or positions[j] < right_lim:
                temp_cs.append(currents[j])

    else: # if 0 < i < len(pos) - 1
        for j in range(len(positions)):
            if positions[j] > left_lim and positions[j] < right_lim:
                temp_cs.append(currents[j])
    
    if len(temp_cs) != 0:
        uniform_currents.append(sum(temp_cs) / len(temp_cs))
        uniform_positions.append(uniform_positions_desired[i])


# Define Fourier variables
L = 2*np.pi                     # Domain  
xs = uniform_positions          # x
fs = uniform_currents           # f(x)

# Compute A0
a0  = 0
for i in range(len(xs) - 1):
    # Trapezoidal integral
    f = (fs[i] + fs[i+1]) / 2
    dt = xs[i+1] - xs[i]
    a0 += f * dt
a0 = a0 / L


# Compute Coefficients
aks, bks = [], []
NUM_COEFFS = 160
for k in range(1, NUM_COEFFS):
    ak = 0
    bk = 0
    for i in range(len(xs) - 1):
        f = (fs[i] + fs[i+1]) / 2
        dt = xs[i+1] - xs[i]
        ak += (2/L) * f*np.cos((2*pi*k*xs[i])/L) * dt
        bk += (2/L) * f*np.sin((2*pi*k*xs[i])/L) * dt
    
    aks.append(ak)
    bks.append(bk)

# Reconstruct the fourier function
fourier = []
NUM_POINTS_DESIRED = 3600*2
fourier_xs = np.linspace(0, 2*np.pi - (2*np.pi/NUM_POINTS_DESIRED), NUM_POINTS_DESIRED)
for x in fourier_xs:
    f = a0
    for k in range(len(aks)):
        f += aks[k]*np.cos((k+1)*2*pi*x/L) + bks[k]*np.sin((k+1)*2*pi*x/L)
    fourier.append(f)

# Plot results
# plt.plot(positions, currents, '.', label = "Data")
plt.plot(xs, fs, label = "Mean and Uniform Data")
plt.plot(fourier_xs, fourier, '.' ,label = "Fourier Function")
plt.title(f"Current fit with uniform sampling in space")
plt.xlabel("Position [rad/s]")
plt.ylabel("Current [A]")
plt.grid()
plt.legend()
plt.show()

iqs_cogging = fourier
# Save data in file
f = open("cogging_tableau.txt", "w")
f.write(f"static const int16_t anti_cogging_tableau[{NUM_POINTS_DESIRED}]" + " = {\n")

for i in range(len(iqs_cogging) - 1):
    f.write(f"{round(iqs_cogging[i] * (1<<16))}, ")

f.write(f"{round(iqs_cogging[-1] * (1<<16))}")
f.write("};")
f.close()
