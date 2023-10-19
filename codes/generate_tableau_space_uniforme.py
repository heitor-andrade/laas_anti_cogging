import matplotlib.pyplot as plt
from utils import get_data
from math import pi
import numpy as np
from scipy.optimize import curve_fit

path = "rolling_data/"
filename = "rolling_data2.txt"

times, positions, velocities, iqs, uqs, supplys = get_data(filename, path)

N_TURNS = 1

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

# # plt.plot(iqs, label="Current")
# plt.plot(times*2*pi/period, velocities, label="Velocity")
# plt.plot(positions, velocities, label="Positions")
# plt.plot(positions, times, '.', label="Positions")
# plt.xlabel("Time [s]")
# plt.ylabel("Velocity [rad/s]")
# plt.legend()
# plt.show()


# Sort by position
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


# # Plot results
# plt.plot(xs, fs, '.', label = "Velocity Data")
# plt.plot(xs, np.fft.ifft(np.fft.fft(fs)), ',' ,label = "Numpy fft")
# plt.title(f"Fit velocity with {NUM_COEFFS} Fourier Coefficients with Padding")
# plt.xlabel("Time [s]")
# plt.ylabel("Velocity [rad/s]")
# plt.grid()
# plt.legend()
# plt.show()

# Get velocities uniform in position

NUM_POINTS = 360*2
uniform_positions = np.linspace(0, 2*pi, NUM_POINTS)
uniform_velocities = []
uniform_times = []

for i in range(len(uniform_positions)):
    temp_vs = []
    temp_ts = []
    left_lim = uniform_positions[i-1]  
    if i != len(uniform_positions) - 1:
        right_lim = uniform_positions[i+1]        
    else:
        right_lim = uniform_positions[1]

    if i == 0 or i == len(uniform_positions) - 1:
        for j in range(len(positions)):
            if positions[j] > left_lim or positions[j] < right_lim:
                temp_vs.append(velocities[j])
                temp_ts.append(times[j])

    else: # if 0 < i < len(pos) - 1
        for j in range(len(positions)):
            if positions[j] > left_lim and positions[j] < right_lim:
                temp_vs.append(velocities[j])
                temp_ts.append(times[j])
    
    uniform_velocities.append(sum(temp_vs) / len(temp_vs))
    uniform_times.append(sum(temp_ts) / len(temp_ts))


times = uniform_times
velocities = uniform_velocities
positions = uniform_positions

# Sort by time
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

# Define Fourier variables
L = period*N_TURNS      # Domain  
xs = times              # x
fs = velocities         # f(x)

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
NUM_COEFFS = 30
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

# Plot results
plt.plot(xs, fs, label = "Velocity Data")
# plt.plot(xs, fourier, '.' ,label = "Fourier Function")
plt.plot(xs, np.fft.ifft(np.fft.fft(fs)), '.' ,label = "Numpy fft")
plt.title(f"Velocity fit with uniform sampling in space")
plt.xlabel("Time [s]")
plt.ylabel("Velocity [rad/s]")
plt.grid()
plt.legend()
plt.show()


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

acc = get_acceleration(times)

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

# Compute acceleration each 0.2 degree
ps = np.linspace(0, 2*pi, 360*5)
times_positions = [get_time(p, m, b) for p in ps]
acs = get_acceleration(times_positions)

# plt.plot(positions, acc, '.', label = "Acceleration")
plt.plot(ps, acs, label = "Acceleration")
plt.title(f"Acceleration (Derivative of fourier velocity)")
plt.xlabel("Position [rad]")
plt.ylabel("Acceleration [rad/(s^2)]")
plt.legend()
plt.show()

# Create Cogging Tableau
acc2cur = 0.000196
iqs_cogging = np.array(acc)*acc2cur*(-1)

plt.plot(positions, iqs_cogging, '.')
plt.title(f"Anti-Cogging Current")
plt.xlabel("Position [rad]")
plt.ylabel("Current [Amperes]")
plt.legend()
plt.show()

# Save data in file
f = open("cogging_tableau.txt", "w")
f.write("float anti_cogging_tableau[1800] = {\n")

for i in range(len(ps) - 1):
    f.write(f"{round(iqs_cogging[i], 5)}, ")

f.write(f"{round(iqs_cogging[-1], 5)}")
f.write("};")
f.close()
