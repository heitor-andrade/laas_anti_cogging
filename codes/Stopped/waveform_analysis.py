import matplotlib.pyplot as plt
from utils import get_data
from math import pi
import numpy as np

# Filter measures of smaller magnitudes for the normal direction
path = "stopped_data/"
filename = "stopped_data5.txt"
data = get_data(filename, path)

times_raw,positions_raw, velocities_raw, iqs_raw, uqs_raw, supplys_raw, errors_pos_raw = data
times,positions, velocities, iqs, uqs, supplys, errors_pos, duty_cicles = [], [], [], [], [], [], [], []

duty_cicles_raw = [uqs_raw[i]/supplys_raw[i] for i in range(len(uqs_raw))]


for i in range(len(positions_raw)):
    pos = positions_raw[i]%(2*pi)
    dc = uqs_raw[i]/supplys_raw[i]

    if pos not in positions:
        positions.append(pos)
        iqs.append(iqs_raw[i])
        uqs.append(uqs_raw[i])
        supplys.append(supplys_raw[i])
        duty_cicles.append(dc)
    else:
        ind = positions.index(pos)
        if abs(dc) > abs(duty_cicles[ind]):
            iqs[ind] = iqs_raw[i]
            uqs[ind] = uqs_raw[i]
            supplys[ind] = supplys_raw[i]
            duty_cicles[ind] = dc


# Filter measures of smaller magnitudes for the reverse direction
filename1 = "stopped_data4_reverse.txt"
times_raw,positions_raw, velocities_raw, iqs_raw, uqs_raw, supplys_raw, errors_pos_raw = get_data(filename1, path)
times_rev,positions_rev, velocities_rev, iqs_rev, uqs_rev, supplys_rev, errors_pos_rev, duty_cicles_rev = [], [], [], [], [], [], [], []

duty_cicles_raw = [uqs_raw[i]/supplys_raw[i] for i in range(len(uqs_raw))]


for i in range(len(positions_raw)):
    pos = positions_raw[i]%(2*pi)
    dc = uqs_raw[i]/supplys_raw[i]

    if pos not in positions_rev:
        positions_rev.append(pos)
        iqs_rev.append(iqs_raw[i])
        uqs_rev.append(uqs_raw[i])
        supplys_rev.append(supplys_raw[i])
        duty_cicles_rev.append(dc)
    else:
        ind = positions_rev.index(pos)
        if abs(dc) > abs(duty_cicles_rev[ind]):
            iqs_rev[ind] = iqs_raw[i]
            uqs_rev[ind] = uqs_raw[i]
            supplys_rev[ind] = supplys_raw[i]
            duty_cicles_rev[ind] = dc


positions += positions_rev
iqs += iqs_rev
uqs += uqs_rev
supplys += supplys_rev
duty_cicles += duty_cicles_rev

plt.plot(np.array(positions)/(2*pi)*360, duty_cicles, ".")
plt.legend()
plt.title(f"Position X Current")
plt.legend()
plt.show()


# Compute dead time, cogging and current waveform and stiction current
dc_dtss         = [0]*len(positions)
dc_cogging      = [0]*len(positions)
iqs_stiction    = [0]*len(positions)
iqs_cogging     = [0]*len(positions)

dc_max  = [-10]*len(positions)
dc_min  = [+10]*len(positions)
iqs_max = [-10]*len(positions)
iqs_min = [+10]*len(positions)

for i in range(len(positions)):

    # Get max and min duty cicles and currents for each position
    for j in range(len(positions)):
        if positions[i] == positions[j]:
            if duty_cicles[j] > dc_max[i]:
                dc_max[i] = duty_cicles[j]
                iqs_max[i] = iqs[j]
            if duty_cicles[j] < dc_min[i]:
                dc_min[i] = duty_cicles[j]
                iqs_min[i] = iqs[j]

    dc_dtss[i]      = (dc_max[i] - dc_min[i]) / 2
    dc_cogging[i]   = (dc_max[i] + dc_min[i]) / 2
    iqs_stiction[i] = (iqs_max[i] - iqs_min[i]) / 2
    iqs_cogging[i]  = (iqs_max[i] + iqs_min[i]) / 2

plt.plot(positions, iqs_max, ".", label="iqs_max")
plt.plot(positions, iqs_min, ",", label="iqs_min")
plt.legend()
plt.title(f"Position X Current")
plt.legend()
plt.show()



dc_dts_max = max(dc_dtss)

dead_times_st = []
for i in range(len(dc_max)):
    if dc_dts_max > dc_max[i] or -dc_dts_max < dc_min[i]:
        dead_times_st.append(dc_dtss[i])

dc_dts = sum(dead_times_st)/len(dead_times_st)

dcs_stiction = []
uqs_stiction = []

# for i in range(len(dc_dtss)):
#     if dc_dts_max < dc_min[i] or (-dc_dts_max) > dc_max[i]:
#         dcs_stiction.append(dc_dtss[i])
#         uqs_stiction.append(dc_dtss[i] * supplys[i])


# uq_stiction = sum(uqs_stiction)/len(uqs_stiction)
# iq_stiction = sum(iqs_stiction)/len(iqs_stiction)
# dc_stiction = sum(dcs_stiction)/len(dcs_stiction)
# dc_dt = dc_dts - dc_stiction


# plt.plot(positions, iqs, ".")
# plt.plot(positions, [iqs_cogging]*len(positions), ".")
# plt.plot(positions, iqs_cogging, ".", label="iqs_cogging")
plt.plot(positions, iqs_stiction, ",", label="iqs_stiction")

# iqs_anti = [iqs_cogging[i] + ]
# plt.plot(positions, iqs_stiction, ",", label="iqs_stiction")


plt.legend()
plt.title(f"Position X Current")
plt.legend()
plt.show()

f = open("anti_cogging.txt", "w")
f.write("Position, I Cogging, I stiction\n")

for i in range(len(times)):
    f.write(str(positions[i]) + ", " 
            + str(iqs_cogging[i]) + ", "
            + str(iqs_cogging[i]) + ", " 
            + "\n")
f.close()
