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


currents_positions, currents_max, currents_min = [], [], []
for i in range(len(positions)):
    temp_currents = []
    for j in range(len(positions_reverse)):
        if abs(positions[i] - positions_reverse[j]) < 2*pi / (5000*2):
            temp_currents.append(currents_reverse[j])

    if len(temp_currents) >= 1:
        temp_currents.append(currents[i])
        currents_positions.append(positions[i])
        currents_max.append(max(temp_currents))
        currents_min.append(min(temp_currents))


# calculate Ist
ists = [(currents_max[i] - currents_min[i])/2 for i in range(len(currents_max))]
ist = sum(ists) / len(ists)

print(ist)

# Plot results
plt.plot(positions + positions_reverse, currents + currents_reverse, '.' ,label = "Actual data")
plt.plot(currents_positions, currents_min, '.' ,label = "Min Currents")
plt.plot(currents_positions, currents_max, '.' ,label = "Max Currents")
plt.plot(currents_positions, [ist]*len(currents_positions), '.' ,label = "Ist current")
plt.title(f"Current X Position using Stopped Algorithm")
plt.xlabel("Position [rad]")
plt.ylabel("Current [A]")
plt.grid()
plt.legend()
plt.show()
