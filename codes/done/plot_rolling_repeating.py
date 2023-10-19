import matplotlib.pyplot as plt
from utils import get_data
from math import pi
import numpy as np
from scipy.optimize import curve_fit

path = "rolling_data/"
filename = "rolling_data2.txt"

times, positions, velocities, iqs, uqs, supplys = get_data(filename, path)

N_TURNS = 3

# Adjust position and time origins
init_time = 0
init_index = 0

for i in range(len(times)):
    if positions[i+1]%(2*pi) == 0:
        init_time = times[i+1]
        init_index = i+1
        break
    elif positions[i]%(2*pi) > 2*pi - 0.1 and positions[i+1]%(2*pi) < 0.1:
        init_time = (times[i] + times[i+1]) / 2
        init_index = i+1
        break

times = times[init_index:]
velocities = velocities[init_index:]
positions = positions[init_index:]

times = [times[i] - init_time for i in range(len(times))]


# Compute time period for one turn 
ds = (positions[-1] - positions[0])         # rad
dt = times[-1] - times[0]                   # s
velocity = ds/dt                            # rad/s
period =  2 * pi/velocity                   # s/(2*pi*rad)

# Union all data velocites in one time 
times = [times[i]%(period*N_TURNS) for i in range(len(times))]
# Limit position range in [0, 2pi]
positions = [positions[i]%(2*pi*N_TURNS) for i in range(len(positions))]

# clean position and time undesired
for i, p in enumerate(positions):
    if p > 5 and times[i] < 0.02:
        positions.pop(i)
        times.pop(i)
        velocities.pop(i)

# Sort time, velocity and position arrays
n = len(times)
for i in range(n):
    change = False
    for j in range(0, n-i-1):
        if times[j] > times[j+1]:
            times[j], times[j+1] = times[j+1], times[j]
            velocities[j], velocities[j+1] = velocities[j+1], velocities[j]
            positions[j], positions[j+1] = positions[j+1], positions[j]
            change = True
    if not change:
        break


# positions = np.array(positions)
# velocities = np.array(velocities)
# times = np.array(times)


# # plt.plot(iqs, label="Current")
# plt.plot(times*2*pi/period, velocities, label="Velocity")
# plt.plot(positions, velocities, label="Positions")
# plt.plot(positions, times, '.', label="Positions")
# plt.xlabel("Time [s]")
# plt.ylabel("Velocity [rad/s]")
# plt.legend()
# plt.show()

# Define Fourier variables
L = period*N_TURNS      # Domain  
xs = times              # x
fs = velocities         # f(x)

NUM_PADDING = 0

# Compute A0
a0  = 0
for i in range(len(xs) - 1):
    # Trapezoidal integral
    f = (fs[i] + fs[i+1]) / 2
    dt = xs[i+1] - xs[i]
    a0 += f * dt
a0 = a0 / L


padding = [a0] * NUM_PADDING
fs = padding + fs + padding
dx = xs[1] - xs[0]
xs = [xs[0] - dx*(i+1) for i in range(NUM_PADDING, 0, -1)] + xs + [xs[-1] + dx*(i+1) for i in range(NUM_PADDING)]

L = L + dx*NUM_PADDING*2

# Compute Coefficients
aks, bks = [], []
NUM_COEFFS = 42
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

# Plot results
plt.plot(xs, fs, label = "Velocity Data")
plt.plot(xs, fourier, label = "Fourier Function")
plt.title(f"Fit velocity with {NUM_COEFFS} Fourier Coefficients")
plt.xlabel("Time [s]")
plt.xlim((period, 2*period))
plt.ylim((50, 51.5))
plt.ylabel("Velocity [rad/s]")
plt.legend()
plt.show()

# plt.plot(positions, acc, '.', label = "Acceleration")
# plt.plot(np.array(xs)*2*pi/period, acc, '.', label = "Acceleration Tempo*2*pi/Periodo")
# # plt.plot(xs, np.gradient(fourier, xs), ',', label = "np.gradient(fourier)")
# plt.title(f"Acceleration (Derivative of fourier velocity)")
# plt.xlabel("Time [s]")
# plt.ylabel("Acceleration [rad/(s^2)]")
# plt.legend()
# plt.show()



# DEFINE RELATION ACCELERATION(POSITION)

# Sort time and velocity array
n = len(positions)
for i in range(n):
    change = False
    for j in range(0, n-i-1):
        if positions[j] > positions[j+1]:
            positions[j], positions[j+1] = positions[j+1], positions[j]
            velocities[j], velocities[j+1] = velocities[j+1], velocities[j]
            times[j], times[j+1] = times[j+1], times[j]
            change = True
    if not change:
        break


# define m, b
def get_time(x, m, b):
    return m*x + b
[m, b], pcov = curve_fit(get_time, positions, times)

# Compute acceleration each 0.1 degree
ps = np.linspace(0, 2*pi, 360*10)
times_positions = [get_time(p, m, b) for p in ps]
acs = get_acceleration(times_positions)

# plt.plot(positions, acc, '.', label = "Acceleration")
# plt.plot(ps, acs, ',', label = "New Acceleration")
# plt.title(f"Acceleration (Derivative of fourier velocity)")
# plt.xlabel("Time [s]")
# plt.ylabel("Acceleration [rad/(s^2)]")
# plt.legend()
# plt.show()










######### USING WHEN TIME AND POSITIONS DO NOT HAVE ADJUST IN THE ORIGIN

# # separate lines and define lim_position
# max_time =0 
# for i in range(len(positions)):
#     if times[i] < 0.04:
#         break
#     max_position_index = i
# lim_position = positions[max_position_index]

# times1          = times[:max_position_index]
# times2          = times[max_position_index+1:]
# positions1      = positions[:max_position_index]
# positions2      = positions[max_position_index+1:]


# [m2, b2], pcov = curve_fit(func, positions2, times2)

# def get_time(position):
#     position = position % (2*pi)
#     if position < lim_position:
#         time = m1*position + b1
#     else:
#         time = m2*position + b2
#     return time

# times_fit = [get_time(p) for p in positions]

# plt.plot(positions, times, '.', label="Data")
# plt.plot(positions, times_fit, ',', label="Fit Function")
# # # plt.plot(positions1, times1, '.', label="Positions 1")
# # # plt.plot(positions2, times2, '.', label="Positions 2")
# plt.xlabel("Position [rad]")
# plt.ylabel("Time [s]")
# plt.legend()
# plt.show()

# # Compute acceleration each 0.1 degree

# ps = np.linspace(0, 2*pi, 360*10)
# times_positions = [get_time(p) for p in ps]
# acs = get_acceleration(times_positions)

# plt.plot(positions, acc, '.', label = "Acceleration")
# # plt.plot(np.array(xs)*2*pi/period, acc, '.', label = "Acceleration Tempo*2*pi/Periodo")
# plt.plot(ps, acs, '.', label = "New Acceleration")
# plt.title(f"Acceleration (Derivative of fourier velocity)")
# plt.xlabel("Time [s]")
# plt.ylabel("Acceleration [rad/(s^2)]")
# plt.legend()
# plt.show()
