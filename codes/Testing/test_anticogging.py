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
count = 0
running = True


while running:
    
    now = time.time()
    count += 1
    if count % 50:
        print(
            "position", round(ud.position0, 2),
            "current", round(ud.current0, 2),
            )

    ud.transfer() 
    t +=dt
    while(time.perf_counter()-t<dt):
        pass

ud.stop() # Terminate