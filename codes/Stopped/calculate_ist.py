import matplotlib.pyplot as plt
from math import pi
import numpy as np

import sys
sys.path.append('../') 
from utils import get_data
from scipy.optimize import curve_fit
from scipy import signal

def get_consolidate_data(filename, path="stopped_data/"):

    times,positions, velocities, iqs, uqs, supplys, errors_pos = get_data(filename, path)
    ps_temp, iqs_temp, dcs_temp = [], [], []
    for i in range(len(positions)):
        # Convert negative positions in positive
        if positions[i] > 0:
            pos = abs(positions[i])%(2*pi)
        else:
            pos = (2*pi)-(abs(positions[i])%(2*pi))

        # Consolidate currents of bigger magnitude for the same position
        dc = uqs[i]/supplys[i]
        if pos not in ps_temp:
            ps_temp.append(pos)
            iqs_temp.append(iqs[i])
            dcs_temp.append(dc)
        else:
            ind = ps_temp.index(pos)
            if abs(iqs[i]) > abs(iqs_temp[ind]):
                iqs_temp[ind] = iqs[i]
            if abs(dc) > abs(dcs_temp[ind]):
                dcs_temp[ind] = dc

    return ps_temp, iqs_temp, dcs_temp


positions, currents, duty_cicles = get_consolidate_data("stopped_data5.txt")
positions_reverse, currents_reverse, duty_cicles_reverse = get_consolidate_data("stopped_data4_reverse.txt")

plt.plot(positions, currents, ".", label = "Forward Data")
plt.plot(positions_reverse, currents_reverse, ".", label = "Reverse Data")
plt.legend()
plt.show()

currents_positions, currents_max, currents_min = [], [], []
dcs_positions, dcs_max, dcs_min = [], [], []

for i in range(len(positions)):
    temp_currents = []
    temp_dcs = []
    for j in range(len(positions_reverse)):
        if abs(positions[i] - positions_reverse[j]) < 2*pi / (5000*2):
            temp_currents.append(currents_reverse[j])
            temp_dcs.append(duty_cicles_reverse[j])

    if len(temp_currents) >= 1 and max(temp_dcs) < 0.2 and min(temp_dcs) > -0.2:
        temp_currents.append(currents[i])
        temp_dcs.append(duty_cicles[i])

        currents_positions.append(positions[i])
        currents_max.append(max(temp_currents))
        currents_min.append(min(temp_currents))
        
        dcs_positions.append(positions[i])
        dcs_max.append(max(temp_dcs))
        dcs_min.append(min(temp_dcs))


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
