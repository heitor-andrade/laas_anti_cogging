import sys
sys.path.append('../')  # Adiciona o diretório pai (pasta_principal) ao PATH
from odri_spi_ftdi import SPIuDriver, PID
import matplotlib.pyplot as plt
import numpy as np
import time
import random
import pygame


# CONSTANTS
dt = 0.001
# Initializing Omodri communication 
ud = SPIuDriver(absolutePositionMode=True, waitForInit=True)
ud.transfer()


# Principal loop 
t = time.perf_counter()
init_time = time.time()
running = True

kp, ki, kd = 0.09, 0, 0.0

pid = PID(kp, ki, kd, 2)
reference = 2*np.pi

times, positions, velocities, velocities_reference, iqs, function_activated = [], [], [], [], [], []
ud.kp0 = 1

while running:
    
    now = time.time() - init_time
    ud.refCurrent0 = pid.compute(ud.velocity0, 0, reference)
    # ud.refCurrent0 =  pid.compute(ud.position0, ud.velocity0, 5.5)

    times.append(round(now, 4))
    positions.append(round(ud.position0, 5))
    velocities.append(round(ud.velocity0, 5))
    velocities_reference.append(round(reference, 5))
    iqs.append(round(ud.current0, 5))
    function_activated.append(ud.kp0)


    if now < 8:
        ud.kp0 = 1
    elif now < 16:
        ud.kp0 = 0
    else:
        running = False

    print(ud.current0, ud.velocity0)
    ud.transfer() 
    t +=dt
    while(time.perf_counter()-t<dt):
        pass

ud.stop() # Terminate

# Save data in file
f = open(f"speed_data/speed_test.txt", "w")
f.write("Time, Position, Velocity, Reference Velocity, Iq, Function Activated\n")

for i in range(len(times)):
    f.write(str(times[i]) + ", " 
            + str(positions[i]) + ", "
            + str(velocities[i]) + ", " 
            + str(velocities_reference[i]) + ", " 
            + str(iqs[i]) + ", "
            + str(function_activated[i])
            + "\n")
f.close()
