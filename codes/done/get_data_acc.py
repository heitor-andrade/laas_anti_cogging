import sys
sys.path.append('../')  # Adiciona o diretório pai (pasta_principal) ao PATH
from odri_spi_ftdi import SPIuDriver, PID
import matplotlib.pyplot as plt
import numpy as np
from math import pi
import time
import random

# CONSTANTS


# Initializing Omodri communication 
ud = SPIuDriver(absolutePositionMode=False, waitForInit=True)
ud.transfer()
times, positions, velocities, iqs, vqs, vsupplys = [], [], [], [], [], []

# Current caracteristics
f = 2       # Hz
A = 0.4    # Amplitude
offset = 0.35

# Principal loop 
dt          = 0.001  # Period Loop
case        = 1
TIME_TEST   = 30
t           = time.perf_counter()
init_time   = time.time()
count = 0

while time.time() - init_time < TIME_TEST:
    
    now = time.time()
    # Set current format
    if case == 1:
        ud.refCurrent0 = A * np.sin(f*2*pi*now)
        
    # Store data while rolling the motor
    if now - init_time > 3:
        times.append(round(now - init_time, 4))
        positions.append(round(ud.position0, 5))
        velocities.append(round(ud.velocity0, 5))
        iqs.append(round(ud.current0, 5))
        vqs.append(round(ud.tension0, 5))
        vsupplys.append(round(ud.supply0, 5))

        print("iq: ", round(ud.current0, 4), 
                "iq reference: ", round(ud.refCurrent0, 4),
                "speed: ", round(ud.velocity0))
        
        # fail safe
        if ud.velocity0 > 160:
            ud.stop()
            print("Stopped!")
            break
        else:
            count = 0

    ud.transfer()
    t += dt
    while(time.perf_counter()-t<dt):
        pass

ud.stop() # Terminate

# Save data in file
f = open("acc_current.txt", "w")
f.write("Time, Position, Velocity, Iq, Uq, Vsupply\n")

for i in range(len(times)):
    f.write(str(times[i]) + ", " 
            + str(positions[i]) + ", "
            + str(velocities[i]) + ", " 
            + str(iqs[i]) + ", "
            + str(vqs[i]) + ", " 
            + str(vsupplys[i])
            + "\n")
f.close()

sys.exit()