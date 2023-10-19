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

filename = "stopped_data_test.txt"
times,positions, velocities, iqs, uqs, supplys, errors_pos = get_data(filename, path)

fig, ax = plt.subplots(2)

index_max = positions.index(max(positions))

positions_reverse = positions[index_max:]
iqs_reverse = iqs[index_max:]

ax[0].plot(positions, iqs, '.', label = "Current Forward")
ax[0].plot(positions_reverse, iqs_reverse, '.', label = "Current Reverse")
# ax[0].plot(velocitiesr, iqsr, '.', label = "Current")
# ax[0].plot(timesr, positions_refr, '.', label = "Reference Position")
# ax[0].plot(times, positions, '.', label = "Actual Position")
ax[0].set_title(f"Time X Position using Stopped Algorithm")
ax[0].set_xlabel("Time [s]")
ax[0].set_ylabel("Position [rad]")
ax[0].grid()
ax[0].legend()

ax[1].plot(timesr, positionsr, '.', label = "Actual position")
ax[1].plot(timesr, positions_refr, '.', label = "Reference position")
# ax[1].plot(positionsr, positions_refr, '.', label = "Velocity")


plt.show()
