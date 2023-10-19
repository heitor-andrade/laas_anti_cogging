import sys
sys.path.append('../')  # Adiciona o diretório pai (pasta_principal) ao PATH
from odri_spi_ftdi import SPIuDriver, PID
import matplotlib.pyplot as plt
import numpy as np
from math import pi
import time
import random

# CONSTANTS
dt = 0.001                                  # Period Loop
UQ_RESOLUTION = 0.000001      # Tension uq resolution


# Initializing Omodri communication 
ud = SPIuDriver(absolutePositionMode=True, waitForInit=True)
ud.transfer()

times, positions, velocities, iqs, vqs, vsupplys = [], [], [], [], [], []

# Principal loop 
t = time.perf_counter()
case = 1
uq_cogmax = 0.045
last_position_stopped = 0
count = 0

while True:
    
    now = time.time()
    # Compute uq minimum needed to surpass dead time, stiction and cogging
    if case == 1:
        if ud.velocity0 == 0:
            last_position_stopped = ud.position0
            uq_cogmax += UQ_RESOLUTION
            print("uq_cogmax:", round(uq_cogmax, 5), "ud.position0:", round(ud.position0, 2))
        elif abs(ud.position0 - last_position_stopped) >= 10*(2*pi):
            case = 2
            uq_cogmin = uq_cogmax - UQ_RESOLUTION
            last_position_moving = ud.position0
            last_time_moving = time.time()

            UQ_RESOLUTION = UQ_RESOLUTION/100
            print("**********GOING TO CASE 2 *******************")
        else:
            print("uq_cogmax defined = ", uq_cogmax)


        ud.refVelocity0 = uq_cogmax

    # Compute uq minimum to maintain the motor rolling
    if case == 2:
        if ud.velocity0 != 0:
            last_time_moving = time.time()
            print("uq_cogmin:", uq_cogmin)
            if abs(ud.position0 - last_position_moving) > 8*pi:
                uq_cogmin -= UQ_RESOLUTION
                last_position_moving = ud.position0

        elif time.time() - last_time_moving > 6:
            case = 3
            init_time = time.time() 
            uq_cogmin += UQ_RESOLUTION
            print("uq_cogmin:", uq_cogmin, "uq_cogmax:", uq_cogmax)
            print("**********GOING TO CASE 3 *******************")


        ud.refVelocity0 = uq_cogmin
    
    # Start rolling the motor
    if case == 3:
        ud.refVelocity0 = uq_cogmax
        print("Rolling applying uq_cogmax(6s). Time = ", round(time.time() - init_time, 1))

        if time.time() - init_time > 2:
            case = 4
            init_time = time.time()
            print("**********GOING TO CASE 4 *******************")

    
    # Roll the motor in the minimum uq possible
    if case == 4:
        ud.refVelocity0 = uq_cogmin
        print("Rolling applying uq_cogmin(6s). Time = ", round(time.time() - init_time, 1))

        if time.time() - init_time > 3:
            case = 5
            num_turns = 30
            init_position = ud.position0
            init_time = time.time()

    # Store data while rolling the motor
    if case == 5:
        ud.refVelocity0 = uq_cogmin
        times.append(round(now - init_time, 4))
        positions.append(round(ud.position0, 5))
        velocities.append(round(ud.velocity0, 5))
        iqs.append(round(ud.current0, 5))
        vqs.append(round(ud.tension0, 5))
        vsupplys.append(round(ud.supply0, 5))

        print(f"Getting data. Turns =  {round(abs(ud.position0 - init_position) / (2*pi))}")
        if abs(ud.position0 - init_position)/(2*pi) > num_turns:
            break

    # fail safe
    if ud.velocity0 > 100: 
        ud.stop()
        print("Stopped!")
        break

    ud.transfer()
    t += dt
    while(time.perf_counter()-t<dt):
        pass

ud.stop() # Terminate

print("uq_cogmax, uq_cogmin, difference", uq_cogmax, uq_cogmin, uq_cogmax - uq_cogmin)

# Save data in file
f = open("rolling_data.txt", "w")
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