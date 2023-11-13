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

# PLOT Velocity

positions = [positions[i] % (2*np.pi) for i in range(len(positions))]

cogging_index = function.index(1)

times_nocog = times[:cogging_index]
velocities_nocog = velocities[:cogging_index]
positions_nocog = positions[:cogging_index]
iqs_nocog = iqs[:cogging_index]

print(len(positions_nocog), len(times_nocog))

velocities_nocog = [velocities_nocog[i] for i in range(len(velocities_nocog)) if times_nocog[i] > 1]
positions_nocog = [positions_nocog[i] for i in range(len(positions_nocog)) if times_nocog[i] > 1]
iqs_nocog = [iqs_nocog[i] for i in range(len(iqs_nocog)) if times_nocog[i] > 1]
times_nocog = [times_nocog[i] for i in range(len(times_nocog)) if times_nocog[i] > 1]

times_cog = np.array(times[cogging_index:]) - times[cogging_index]
velocities_cog = velocities[cogging_index:]
positions_cog = positions[cogging_index:]
iqs_cog = iqs[cogging_index:]

velocities_cog = [velocities_cog[i] for i in range(len(velocities_cog)) if times_cog[i] > 1]
positions_cog = [positions_cog[i] for i in range(len(positions_cog)) if times_cog[i] > 1]
iqs_cog = [iqs_cog[i] for i in range(len(iqs_cog)) if times_cog[i] > 1]
times_cog = [times_cog[i] for i in range(len(times_cog)) if times_cog[i] > 1]


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

print("Mean squared error without and with anticogging: ", round(mse1, 2), round(mse2, 2))



