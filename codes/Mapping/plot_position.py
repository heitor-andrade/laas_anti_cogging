import matplotlib.pyplot as plt
from math import pi
import numpy as np

import sys
sys.path.append('../') 
from utils import get_data
from scipy.optimize import curve_fit
from scipy import signal

filename = "stopped_data_ref.txt"
path = "./"
timesr, positionsr, positions_refr, velocitiesr, iqsr, uqs, supplys = get_data(filename, path)

filename = "stopped_data.txt"
times,positions, velocities, iqs, uqs, supplys, errors_pos = get_data(filename, path)

fig, ax = plt.subplots(2)

index_max = positions.index(max(positions))

positions_reverse = positions[index_max:]
iqs_reverse = iqs[index_max:]

ax[0].plot(positions, iqs, '.', label = "Current Forward")
ax[0].plot(positions_reverse, iqs_reverse, '.', label = "Current Reverse")
ax[0].set_title(f"Position X Current")
ax[0].set_xlabel("Position [rad]")
ax[0].set_ylabel("Current [A]")
ax[0].grid()
ax[0].legend()

ax[1].plot(np.array(timesr)/60, positionsr, label = "Actual position")
ax[1].plot(np.array(timesr)/60, positions_refr, ':', label = "Reference position")
ax[1].set_title(f"Time X Position")
ax[1].set_xlabel("Time [min]")
ax[1].set_ylabel("Position [rad]")
ax[1].grid()
ax[1].legend()

plt.show()
