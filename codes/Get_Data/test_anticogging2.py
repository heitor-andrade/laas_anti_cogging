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

# PYGAME
pygame.init()
WIDTH, HEIGHT = 600, 400
screen = pygame.display.set_mode((WIDTH, HEIGHT))
pygame.display.set_caption("Simulação de Temperatura")
WHITE = (255, 255, 255)
BLACK = (0, 0, 0)
FONT_COLOR = BLACK


kp, ki, kd = 0.09, 0, 0.0

pid = PID(kp, ki, kd, 2)
reference = 2*np.pi

times, positions, velocities, velocities_reference, iqs, function_activated = [], [], [], [], [], []

speed_control = False

while running:
    
    now = time.time() - init_time

    events = pygame.event.get()
    for event in events:
        if event.type == pygame.KEYDOWN:
            if event.key == pygame.K_q:  # Pressionar a tecla "s" para finalizar
                running = False
            elif event.key == pygame.K_a:
                ud.kp0 = 0
                print("Anti-cogging activated !")
            elif event.key == pygame.K_s:
                ud.kp0 = 1
                print("Nothing activated !")
            elif event.key == pygame.K_d:
                # ud.kp0 = 2
                # print("Anti-cogging and Anti Stiction activated !")
                speed_control = True
                init_control = time.time()


    pygame.display.flip()
    screen.fill(WHITE)



    ud.transfer() 
    t +=dt
    while(time.perf_counter()-t<dt):
        pass

ud.stop() # Terminate