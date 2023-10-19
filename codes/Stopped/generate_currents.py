import matplotlib.pyplot as plt
from math import pi
import numpy as np

import sys
sys.path.append('../') 
from utils import get_data
from scipy.optimize import curve_fit
from scipy import signal


path = "stopped_data/"
filenames = ["stopped_data4_reverse.txt", "stopped_data5.txt"]

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

# Sort by position
# n = len(ps)
# for i in range(n):
#     change = False
#     for j in range(0, n-i-1):
#         if ps[j] > ps[j+1]:
#             ps[j], ps[j+1] = ps[j+1], ps[j]
#             cs[j], cs[j+1] = cs[j+1], cs[j]
#             dcs[j], dcs[j+1] = dcs[j+1], dcs[j]
#             change = True
#     if not change:
#         break

# Combine os arrays em pares (posição, corrente)
combined = list(zip(ps, cs))
sorted_combined = sorted(combined, key=lambda x: x[0])
positions, currents = zip(*sorted_combined)
positions = list(positions)
currents = list(currents)


# Separate the intervals

# Clean currents
MINIMUM_INTERVAL = 0.015         # rad
changed = False
signal_currents = -1         # 1 for positive and -1 for negatives
i = 0
while i != (len(positions) - 1):

    if currents[i] < 0:
        signal_currents = -1
    else:
        signal_currents = 1

    if currents[i+1] < -0.2 and signal_currents == 1:
        changed = True

    if changed:
        changed = False
    else:
        if abs(currents[i+1] - currents[i]) > 0.05:
            currents.pop(i+1)
            positions.pop(i+1)

    if positions[i+1] - positions[i] > MINIMUM_INTERVAL:
        if positions[i] - positions[i-1] > MINIMUM_INTERVAL:
            positions.pop(i)
            currents.pop(i)

    i+=1

i = 0
while i != (len(positions) - 1):
    if abs(currents[i+1] - currents[i]) > 0.05:
        currents.pop(i+1)
        positions.pop(i+1)
    i+=1

# Get end intervals
positions_end = []
currents_end = []
index_end = []
for i in range(len(positions)-1):
    if positions[i+1] - positions[i] > MINIMUM_INTERVAL:
        positions_end.append(positions[i])
        currents_end.append(currents[i])
        index_end.append(i)

positions_end.append(positions[-1])
currents_end.append(currents[-1])

# positions = positions[:index_end[0]]
# currents = currents[:index_end[0]]
def cubic(x, a, b, c, d):
    return a*x**3 + b*x**2 + c*x + d
def line(x, m, n):
    return m*x + n

fpositions, fcurrents = [], []
uniform_currents = []

dx = 0.001

for i in range(len(index_end)):
    # Get interval
    if i == 0:
        xs = positions[0: index_end[0] + 1]
        ys = currents[0: index_end[0] + 1]
    else:
        xs = positions[index_end[i-1] + 1: index_end[i]]
        ys = currents[index_end[i-1] + 1: index_end[i]]

    # Interpolate Cubic
    [a, b, c, d], _ = curve_fit(cubic, xs, ys)
    xs = np.arange(xs[0], xs[-1] + dx, dx)
    ys = [cubic(x, a, b, c, d) for x in xs]

    fpositions += list(xs)
    fcurrents += list(ys)

    # Get space between intervals
    xs = positions[index_end[i] - 2 : index_end[i] + 4]
    ys = currents[index_end[i] - 2 : index_end[i] + 4]

    # Interpolate line
    [m, n], _ = curve_fit(line, xs, ys)
    xs = np.arange(fpositions[-1] + dx, xs[-3] + dx, dx)
    ys = [line(i, m, n) for i in xs]

    # Windowing
    N = len(xs)
    tukey = signal.windows.tukey(N, alpha=0.8)
    ys = ys * (1 - tukey)

    fpositions += list(xs)
    fcurrents += list(ys)

# Get last interval
xs = positions[index_end[-1] + 1:]
ys = currents[index_end[-1] + 1:]

# Interpolate Cubic
[a, b, c, d], _ = curve_fit(cubic, xs, ys)
xs = np.arange(xs[0], xs[-1] + dx, dx)
ys = [cubic(x, a, b, c, d) for x in xs]

fpositions += list(xs)
fcurrents += list(ys)

# Link last and first interval
delta_y = fcurrents[0] - fcurrents[-1]
delta_x = (2*pi - fpositions[-1]) + (fpositions[0])
dy = delta_y / delta_x

pos = fpositions[-1]
changed_side = False


while pos <= 2*pi:
    pos += dx
    fcurrents.append(fcurrents[-1] + (dx * dy))
    fpositions.append(pos)

pos = fpositions[0]
while pos >= 0:
    pos -= dx
    fcurrents.insert(0, fcurrents[0] + (-dx * dy))
    fpositions.insert(0, pos)


# Get Uniform positions and currents
uniform_positions = np.arange(0, 2*pi + dx, dx)
# Sample the current uniformly in space
uniform_currents = []
count = 0

for i in range(len(uniform_positions)):
    temp_currents = []

    left_lim = uniform_positions[i-1]  
    if i != len(uniform_positions) - 1:
        right_lim = uniform_positions[i+1]        
    else:
        right_lim = uniform_positions[1]

    if i == 0 or i == len(uniform_positions) - 1:
        for k in range(len(fpositions)):
            if fpositions[k] > left_lim or fpositions[k] < right_lim:
                temp_currents.append(fcurrents[k])

    else:
        for k in range(len(fpositions)):
            if fpositions[k] > left_lim and fpositions[k] < right_lim:
                temp_currents.append(fcurrents[k])
    
    if len(temp_currents) != 0:
        uniform_currents.append(sum(temp_currents) / len(temp_currents))
    else:
        if len(uniform_currents) != 0:
            uniform_currents.append(uniform_currents[-1])
        else:
            uniform_currents.append(0)
        count +=1


# Plot results
plt.plot(positions, currents,'.' ,label = "Actual data")
plt.plot(np.array(uniform_positions)*360/(2*pi), uniform_currents, '.' ,label = "Uniform Data")
# plt.plot(fpositions, fcurrents, '.' ,label = "Final Function")
# plt.title(f"Current fit with uniform sampling in space")
plt.title(f"Current X Position using Stopped Algorithm")
plt.xlabel("Position [rad]")
plt.ylabel("Current [A]")
plt.grid()
plt.legend()
plt.show()




# Save data in file

NUM_POINTS = len(uniform_positions)
f = open("cogging_tableau_CF.txt", "w")
f.write(f"static const int16_t anti_cogging_tableau[{NUM_POINTS}]" + "= {\n")

for i in range(len(uniform_positions) - 1):
    u = int(uniform_currents[i] * (1<<16))
    f.write(f"{round(u)}, ")

u = int(uniform_currents[-1] * (1<<16))
f.write(f"{round(u)}")
f.write("};")
f.close()
