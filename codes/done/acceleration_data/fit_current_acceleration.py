import matplotlib.pyplot as plt
from math import pi
import numpy as np
from scipy.optimize import curve_fit

import sys
sys.path.append('../') 
from utils import get_data

path = "./"
filename = "acc_current.txt"

times, positions, velocities, iqs, uqs, supplys = get_data(filename, path)

# Take out the initial 5 seconds
time_init = 5
positions = [positions[i] for i in range(len(times)) if times[i] > time_init]
velocities = [velocities[i] for i in range(len(times)) if times[i] > time_init]
iqs = [iqs[i] for i in range(len(times)) if times[i] > time_init]
times = [times[i] for i in range(len(times)) if times[i] > time_init]

# Fit Velocity

# plt.plot(times, velocities, label="data")
# plt.plot(times, fit, label="fit")
# plt.legend()
# plt.show()


# Define Fourier variables
L = times[-1] - times[0]        # Domain  
xs = times                      # x
fs = velocities                 # f(x)


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
NUM_COEFFS = 60
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

# Reconstruct the function velocity
fourier = []
for i in range(len(xs)):
    f = a0
    x = xs[i]
    for k in range(len(aks)):
        f += aks[k]*np.cos((k+1)*2*pi*x/L) + bks[k]*np.sin((k+1)*2*pi*x/L)
    fourier.append(f)


# # Plot results
# plt.plot(xs, fs, label = "Velocity Data")
# plt.plot(xs, fourier, ',', label = "Fourier Function")
# plt.title(f"Fit velocity with {NUM_COEFFS} Fourier Coefficients")
# plt.xlabel("Time [s]")
# plt.ylabel("Velocity [rad/s]")
# plt.legend()
# plt.show()


# Compute acceleration = derivative of velocity
def get_acceleration(times):
    acc = []
    for t in times:
        a = 0
        for k in range(len(aks)):
            # f += aks[k]*np.cos((k+1)*2*pi*x/L) + bks[k]*np.sin((k+1)*2*pi*x/L)
            a += -aks[k]* ((k+1)*2*pi/L)*np.sin((k+1)*2*pi*t/L) + bks[k]* ((k+1)*2*pi/L) * np.cos((k+1)*2*pi*t/L)
        acc.append(a)
    return acc

acc = get_acceleration(xs)

# plt.plot(times, acc, '.', label = "Acceleration")
# plt.plot(times, np.gradient(fourier, times), ',', label = "np.gradient(fourier)")
# plt.title(f"Acceleration (Derivative of fourier velocity)")
# plt.xlabel("Time [s]")
# plt.ylabel("Acceleration [rad/(s^2)]")
# plt.legend()
# plt.show()


# Least squares current x acceleration
def func(x, a):
    return a*x

# acc_scale = max(iqs)/max(acc)

[acc_scale], _ = curve_fit(func, acc, iqs)

print(acc_scale)

iqs_fit = [func(a, acc_scale) for a in acc]

# plt.plot(iqs, acc, '.')
# plt.plot(times, np.array(acc)*acc_scale, '.', label="acceleration")
plt.plot(times, iqs, '.', label="Data")
plt.plot(times, iqs_fit, '-', label=f"Current = {round(acc_scale, 6)}*Acceleration")
plt.title(f"Fit acceleration to the current")
plt.xlabel("Time [s]")
plt.ylabel("Current [A]")
plt.legend()
plt.show()

