import matplotlib.pyplot as plt
from math import pi
import numpy as np

import sys
sys.path.append('../') 
from utils import get_data

path = "stopped_data/"
filenames = ["stopped_data4_reverse.txt", "stopped_data5.txt"]
# filenames = [ "stopped_data5.txt"]

DEG2RAD = 2*pi/360

ps = []
cs = []
dcs = []

for filename in filenames:
    datas = get_data(filename, path)
    if len(datas) == 6:
        times,positions, velocities, iqs, uqs, supplys = datas
    else:
        times,positions, velocities, iqs, uqs, supplys, errors_pos = datas

    positions_degree, iqs_degree, duty_cicles = [], [], []
    for i in range(len(positions)):
        if positions[i] > 0:
            pos = abs(positions[i])%(2*pi)
        else:
            pos = (2*pi)-(abs(positions[i])%(2*pi))

        pos = pos
        dc = uqs[i]/supplys[i]
        if pos not in positions_degree:
            positions_degree.append(pos)
            iqs_degree.append(iqs[i])
            duty_cicles.append(dc)
        else:
            ind = positions_degree.index(pos)
            if abs(iqs[i]) > abs(iqs_degree[ind]):
                iqs_degree[ind] = iqs[i]
            if abs(dc) > abs(duty_cicles[ind]):
                duty_cicles[ind] = iqs[i]

    ps += positions_degree
    cs += iqs_degree
    dcs += duty_cicles

    # plt.plot(positions_degree, duty_cicles, ".", label = filename)
    # plt.plot(np.array(positions_degree), iqs_degree, ".", label = "current")
    # plt.plot(np.array(positions_degree), [duty_cicles[i] - iqs_degree[i] for i in range(len(duty_cicles))], ".", label = "current")
    # plt.plot(duty_cicles, iqs_degree, ".", label = "dcxcurrent")


# # plt.plot(ps, dcs, ".", label = "Duty cicle")
# plt.plot(ps, cs, ".", label = "Current")
# plt.title(f"Position X Current")
# plt.legend()
# plt.show()



# Sample the current uniformly in space
NUM_POINTS = 360*20
uniform_positions_init = np.linspace(0, 2*pi, NUM_POINTS)
uniform_positions = []
uniform_currents = []
uniform_dutycicles = []

for i in range(len(uniform_positions_init)):
    temp_cs = []
    temp_dcs = []

    left_lim = uniform_positions_init[i-1]  
    if i != len(uniform_positions_init) - 1:
        right_lim = uniform_positions_init[i+1]        
    else:
        right_lim = uniform_positions_init[1]

    if i == 0 or i == len(uniform_positions_init) - 1:
        for k in range(len(ps)):
            if ps[k] > left_lim or ps[k] < right_lim:
                temp_cs.append(cs[k])
                temp_dcs.append(dcs[k])

    else:
        for k in range(len(ps)):
            if ps[k] > left_lim and ps[k] < right_lim:
                temp_cs.append(cs[k])
                temp_dcs.append(dcs[k])
    
    uniform_positions.append(uniform_positions_init[i])    
    if len(temp_cs) != 0:
        uniform_currents.append(sum(temp_cs) / len(temp_cs))
        uniform_dutycicles.append(sum(temp_dcs) / len(temp_dcs))
    else:
        uniform_currents.append(0)
        uniform_dutycicles.append(0)




def Fourier(L, xs, fs, NUM_COEFFS = 80, num_fourier_points = None):
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
    for k in range(1, NUM_COEFFS):
        ak = 0
        bk = 0
        dt = (xs[-1] - xs[0])/len(xs)
        for i in range(len(xs) - 1):
            f = (fs[i])
            # dt = xs[i+1] - xs[i]
            ak += (2/L) * f*np.cos((2*pi*k*xs[i])/L) * dt
            bk += (2/L) * f*np.sin((2*pi*k*xs[i])/L) * dt
        
        aks.append(ak)
        bks.append(bk)

    # Reconstruct the fourier function
    if num_fourier_points != None:
        xs = np.linspace(xs[0], xs[-1], num_fourier_points)

    fourier = []
    for i in range(len(xs)):
        f = a0
        x = xs[i]
        for k in range(len(aks)):
            f += aks[k]*np.cos((k+1)*2*pi*x/L) + bks[k]*np.sin((k+1)*2*pi*x/L)
        fourier.append(f)
    
    return fourier

# Define Fourier variables
xs = uniform_positions          # x
fs = uniform_currents           # f(x)
fourier = Fourier(xs[-1] - xs[0], uniform_positions, uniform_currents)

# Plot results
plt.plot(ps, cs,'.' ,label = "Actual data")
# plt.plot(xs, fourier ,label = "Fourier Function")
# plt.plot(uniform_positions, uniform_currents, '.',label = "Uniform Data")
# plt.plot(ifft_result, '.' ,label = "Fourier Function")
# plt.plot(np.fft.ifft(fft_result), '.' ,label = "Numpy fft")
plt.title(f"Current fit with uniform sampling in space")
plt.xlabel("Position [rad]")
plt.ylabel("Current [A]")
plt.grid()
plt.legend()
plt.show()

# iqs_cogging = np.array(fourier)
iqs_cogging = np.array(uniform_currents)

# Save data in file
f = open("cogging_tableau.txt", "w")
f.write(f"static const int16_t anti_cogging_tableau[{NUM_POINTS}]" + "= {\n")

for i in range(len(uniform_positions) - 1):
    u = int(uniform_currents[i] * (1<<16))
    f.write(f"{round(u)}, ")

u = int(uniform_currents[-1] * (1<<16))
f.write(f"{round(u)}")
f.write("};")
f.close()
