import matplotlib.pyplot as plt
from math import pi
import numpy as np

import sys
sys.path.append('../') 
from utils import get_data
from scipy.optimize import curve_fit
from scipy import signal

filename = "speed_test.txt"
path = "./"
times, positions, velocities, velocities_reference, iqs, function = get_data(filename, path)

# PLOT Velocity

positions = [positions[i] % (2*np.pi) for i in range(len(positions))]
cogging_index = function.index(0)
times_nocog = times[:cogging_index]
velocities_nocog = velocities[:cogging_index]
positions_nocog = positions[:cogging_index]
iqs_nocog = iqs[:cogging_index]

times_cog = np.array(times[cogging_index:]) - times[cogging_index]
velocities_cog = velocities[cogging_index:]
positions_cog = positions[cogging_index:]
iqs_cog = iqs[cogging_index:]

# plt.plot(times_nocog, iqs_nocog, '.', label = "Without Anti-Cogging")
# plt.plot(times_cog, iqs_cog, '.', label = "With Anti-Cogging")
# plt.title(f"Time X Current, kp = 0.09")
# plt.xlabel("Time [s]")
# plt.ylabel("Current [A]")
# plt.grid()
# plt.legend()
# plt.show()

fig, ax = plt.subplots(2)
ax[0].plot(times_nocog, velocities_nocog, '.', label = "Without Anti-Cogging")
ax[0].plot(times_cog, velocities_cog, '.', label = "With Anti-Cogging")
ax[0].plot(times_cog, len(times_cog) * [velocities_reference[0]], label = "Reference")
# ax[0].plot(times, velocities_reference, '.', label = "Reference")
ax[0].set_title(f"Speed control with kp = 0.09")
ax[0].set_xlabel("Time [s]")
ax[0].set_ylabel("Velocity [rad/s]")
ax[0].grid()
ax[0].legend()

error_nocog = [velocities_nocog[i] - velocities_reference[i] for i in range(len(velocities_nocog))]
error_cog   = [velocities_cog[i] - velocities_reference[i] for i in range(len(velocities_cog))]
ax[1].plot(positions_nocog, error_nocog , '.', label = "Without Anti-Cogging")
ax[1].plot(positions_cog, error_cog, '.', label = "With Anti-Cogging")
ax[1].set_title(f"Speed error over position")
ax[1].set_xlabel("Position [rad/s]")
ax[1].set_ylabel("Error Velocity [rad/s]")
ax[1].legend()
ax[1].grid()
plt.show()

mse1 = np.mean(np.array(error_nocog)**2)
mse2 = np.mean(np.array(error_cog)**2)

print(mse1, mse2)

# PLOT POSITION

# fig, ax = plt.subplots(2)
# cogging_index = function.index(1)
# times_nocog = times[:cogging_index]
# velocities_nocog = velocities[:cogging_index]
# positions_nocog = np.array(positions[:cogging_index]) - positions[0]

# times_cog = np.array(times[cogging_index:]) - times[cogging_index]
# velocities_cog = velocities[cogging_index:]
# positions_cog = np.array(positions[cogging_index:]) - positions[cogging_index]

# ax[0].plot(times_nocog, positions_nocog, label = "Without Anti-Cogging")
# ax[0].plot(times_cog, positions_cog, label = "With Anti-Cogging")
# ax[0].plot(times_cog, times_cog * 2*np.pi, label = "reference")
# ax[0].set_title(f"Time X Position using Stopped Algorithm")
# ax[0].set_xlabel("Time [s]")
# ax[0].set_ylabel("Position [rad]")
# ax[0].grid()
# ax[0].legend()

# plt.show()

