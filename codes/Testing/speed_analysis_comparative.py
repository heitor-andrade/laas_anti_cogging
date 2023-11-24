import matplotlib.pyplot as plt
from math import pi
import numpy as np

import sys
sys.path.append('../') 
from utils import get_data
from scipy.optimize import curve_fit
from scipy import signal

filename = "speed_test1.txt"
path = "./speed_data/"
times, positions, velocities, velocities_reference, iqs, function = get_data(filename, path)

positions = [positions[i] % (2*np.pi) for i in range(len(positions))]
cogging_index = function.index(1)
times_cog1 = np.array(times[cogging_index:]) - times[cogging_index]
velocities_cog1 = velocities[cogging_index:]
positions_cog1 = positions[cogging_index:]
iqs_cog1 = iqs[cogging_index:]

velocities_cog1 = [velocities_cog1[i] for i in range(len(velocities_cog1)) if times_cog1[i] > 1 and times_cog1[i] < 8]
positions_cog1 = [positions_cog1[i] for i in range(len(positions_cog1)) if times_cog1[i] > 1 and times_cog1[i] < 8]
iqs_cog1 = [iqs_cog1[i] for i in range(len(iqs_cog1)) if times_cog1[i] > 1 and times_cog1[i] < 8]
times_cog1 = [times_cog1[i] for i in range(len(times_cog1)) if times_cog1[i] > 1 and times_cog1[i] < 8]


filename = "speed_test.txt"
times, positions, velocities, velocities_reference, iqs, function = get_data(filename, path)

positions = [positions[i] % (2*np.pi) for i in range(len(positions))]
cogging_index = function.index(0)
times_cog2 = np.array(times[cogging_index:]) - times[cogging_index]
velocities_cog2 = velocities[cogging_index:]
positions_cog2 = positions[cogging_index:]
iqs_cog2 = iqs[cogging_index:]

velocities_cog2 = [velocities_cog2[i] for i in range(len(velocities_cog2)) if times_cog2[i] > 1]
positions_cog2 = [positions_cog2[i] for i in range(len(positions_cog2)) if times_cog2[i] > 1]
iqs_cog2 = [iqs_cog2[i] for i in range(len(iqs_cog2)) if times_cog2[i] > 1]
times_cog2 = [times_cog2[i] for i in range(len(times_cog2)) if times_cog2[i] > 1]


fig, ax = plt.subplots(2)
ax[0].plot(times_cog1, velocities_cog1, '.', label = "Old Anti-Cogging")
ax[0].plot(times_cog2, velocities_cog2, '.', label = "New Anti-Cogging")
ax[0].plot(times_cog1, len(times_cog1) * [velocities_reference[0]], label = "Reference")
# ax[0].plot(times, velocities_reference, '.', label = "Reference")
ax[0].set_title(f"Speed control with kp = 0.09")
ax[0].set_xlabel("Time [s]")
ax[0].set_ylabel("Velocity [rad/s]")
ax[0].grid()
ax[0].legend()

error_cog2 = [velocities_cog2[i] - velocities_reference[i] for i in range(len(velocities_cog2))]
error_cog1   = [velocities_cog1[i] - velocities_reference[i] for i in range(len(velocities_cog1))]
ax[1].plot(positions_cog1, error_cog1, '.', label = "Old Anti-Cogging")
ax[1].plot(positions_cog2, error_cog2 , '.', label = "New Anti-Cogging")
ax[1].set_title(f"Speed error over position")
ax[1].set_xlabel("Position [rad/s]")
ax[1].set_ylabel("Error Velocity [rad/s]")
ax[1].legend()
ax[1].grid()
plt.show()


# fig, ax = plt.subplots(2)
# ax[0].plot(times_nocog, velocities_nocog, '.', label = "Without Anti-Cogging")
# ax[0].plot(times_cog, velocities_cog, '.', label = "With Anti-Cogging")
# ax[0].plot(times_cog, len(times_cog) * [velocities_reference[0]], label = "Reference")
# # ax[0].plot(times, velocities_reference, '.', label = "Reference")
# ax[0].set_title(f"Speed control with kp = 0.09")
# ax[0].set_xlabel("Time [s]")
# ax[0].set_ylabel("Velocity [rad/s]")
# ax[0].grid()
# ax[0].legend()

# error_nocog = [velocities_nocog[i] - velocities_reference[i] for i in range(len(velocities_nocog))]
# error_cog   = [velocities_cog[i] - velocities_reference[i] for i in range(len(velocities_cog))]
# ax[1].plot(positions_nocog, error_nocog , '.', label = "Without Anti-Cogging")
# ax[1].plot(positions_cog, error_cog, '.', label = "With Anti-Cogging")
# ax[1].set_title(f"Speed error over position")
# ax[1].set_xlabel("Position [rad/s]")
# ax[1].set_ylabel("Error Velocity [rad/s]")
# ax[1].legend()
# ax[1].grid()
# plt.show()
