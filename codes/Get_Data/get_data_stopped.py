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
ENCODER_RESOLUTION = 5000 * 4
ENCODER_RESOLUTION_RAD = 2 * np.pi / ENCODER_RESOLUTION

times, positions, velocities, iqs, vqs, vsupplys, errors_pos = [], [], [], [], [], [], []
timesr, positionsr, positions_refr, velocitiesr, iqsr, vqsr, vsupplysr = [], [], [], [], [], [], []


# Pygame INITIALIZATION
pygame.init()
WIDTH, HEIGHT = 600, 400
screen = pygame.display.set_mode((WIDTH, HEIGHT))
pygame.display.set_caption("Anticogging Mapping")
WHITE = (255, 255, 255)
BLACK = (0, 0, 0)
FONT_COLOR = BLACK
font = pygame.font.Font(None, 36)
def render_text(text, x, y):
    text_surface = font.render(text, True, FONT_COLOR)
    text_rect = text_surface.get_rect()
    text_rect.center = (x, y)
    screen.blit(text_surface, text_rect)
key_to_number = {
    pygame.K_KP0: '0',
    pygame.K_KP1: '1',
    pygame.K_KP2: '2',
    pygame.K_KP3: '3',
    pygame.K_KP4: '4',
    pygame.K_KP5: '5',
    pygame.K_KP6: '6',
    pygame.K_KP7: '7',
    pygame.K_KP8: '8',
    pygame.K_KP9: '9',
    pygame.K_PERIOD: '.',  # Ponto decimal
}


# Initializing Omodri communication 
ud = SPIuDriver(absolutePositionMode=True, waitForInit=True)
ud.transfer()


# Principal loop 
t = time.perf_counter()
init_time = time.time()
INTERVAL_POSITION = 0.002
position_ref = -INTERVAL_POSITION
# position_ref = 1*2*np.pi /360


# Position control PID
Kp = 3
Kd = 0.006
Ki = 1

# After initialization: Kp = 30, Ki = ?, Kd = 0.006

pid0 = PID(Kp,Ki,Kd, 0.5)

last_time_moving = 0
stopped = False
direction = 1
count_screen = 0

simulating = True
text = ''
save = False
count_store = 0

while (ud.position0  > -0.01 or time.time() - init_time < 5) and simulating :
    now = time.time()
    error = position_ref - ud.position0


    if ud.velocity0 > 100: # fail safe
        ud.stop()
        print("stopei")
        break

    
    # ud.refCurrent0 = pid0.compute(ud.position0, ud.velocity0, position_ref,0.0)
    ud.kp0 = pid0.Kp
    ud.kd0 = pid0.Kd
    ud.refCurrent0 = pid0.Ki
    ud.refPosition0 = position_ref

    # Store data
    positions_refr.append(round(position_ref, 5))
    positionsr.append(round(ud.position0, 5))
    timesr.append(round(now - init_time, 4))
    velocitiesr.append(round(ud.velocity0, 5))
    iqsr.append(round(ud.current0, 5))
    vqsr.append(round(ud.tension0, 5))
    vsupplysr.append(round(ud.supply0, 5))

    if abs(ud.velocity0) == 0.00 and now - init_time > 5:
        if  stopped == False and abs(ud.position0 - position_ref) < INTERVAL_POSITION: # and now - last_time_moving > dt*3:
            stopped = True
    else:
        stopped = False
        last_time_moving = now

    if(stopped):
        times.append(round(now - init_time, 4))
        positions.append(round(ud.position0, 5))
        velocities.append(round(ud.velocity0, 5))
        iqs.append(round(ud.current0, 5))
        vqs.append(round(ud.tension0, 5))
        vsupplys.append(round(ud.supply0, 5))
        errors_pos.append(round(error, 5))

        stopped = False
        count_store += 1
        if count_store == 5:
            position_ref +=  direction * INTERVAL_POSITION
            count_store = 0

    if abs(position_ref) > 2*np.pi and now > 5:
        direction = -1
    
    ud.transfer() # transfer
    
    ### PYGAME
    # Process input keyboard
    events = pygame.event.get()
    for event in events:
        if event.type == pygame.QUIT:
            simulating = False
        elif event.type == pygame.KEYDOWN:
            if event.key == pygame.K_s:  # Pressionar a tecla "s" para finalizar
                simulating = False
            if event.key in key_to_number:
                text += key_to_number[event.key]
            if event.key == pygame.K_p:
                if text == '':
                    pid0.Kp += 1
                else:
                    pid0.Kp = float(text)
                text = ''
            if event.key == pygame.K_i:
                if text == '':
                    pid0.Ki += 0.01
                else:
                    pid0.Ki = float(text)
                text = ''
            if event.key == pygame.K_d:
                if text == '':
                    pid0.Kd += 0.0001
                else:
                    pid0.Kd = float(text)
                text = ''

    # Renderize screen
    count_screen += 1
    if count_screen == 100:
        screen.fill(WHITE)
        render_text(f"Ref Pos: {position_ref:.3f} rad", WIDTH // 2, HEIGHT // 2 -120)
        render_text(f"Pos: {(ud.position0):.3f} rad", WIDTH // 2, (HEIGHT // 2)-90)
        render_text(f"Error: {(error):.4f} rad", WIDTH // 2, (HEIGHT // 2)-60)
        render_text(f"Current: {(ud.current0):.2f} A", WIDTH // 2, (HEIGHT // 2)-30)
        render_text(f"Input: {text}", WIDTH // 2, HEIGHT // 2)
        render_text(f"Kp: {pid0.Kp}", WIDTH // 2, HEIGHT // 2 + 30)
        render_text(f"Ki: {pid0.Ki}", WIDTH // 2, HEIGHT // 2 + 60)
        render_text(f"Kd: {pid0.Kd}", WIDTH // 2, HEIGHT // 2 + 90)
        render_text(f"Time: {round((now - init_time)/60, 1)} min", WIDTH // 2, HEIGHT // 2 + 120)
        if position_ref != 0:
            render_text(f"Time estimated: {round((round((now - init_time)/60) * 4*np.pi / position_ref)/60, 1) } h", WIDTH // 2, HEIGHT // 2 + 150)
        pygame.display.flip()
        count_screen = 0


    t +=dt
    while(time.perf_counter()-t<dt):
        pass

ud.stop() # Terminate



# Save data in file
f = open("stopped_data_test.txt", "w")
f.write("Time, Position, Velocity, Iq, Uq, Vsupply, Error Position\n")

for i in range(len(times)):
    f.write(str(times[i]) + ", " 
            + str(positions[i]) + ", "
            + str(velocities[i]) + ", " 
            + str(iqs[i]) + ", "
            + str(vqs[i]) + ", " 
            + str(vsupplys[i]) + ", " 
            + str(errors_pos[i])
            + "\n")
f.close()


f = open("stopped_data_ref.txt", "w")
f.write("Time, Position, Reference Position, Velocity, Iq, Uq, Vsupply\n")

for i in range(len(timesr)):
    f.write(str(timesr[i]) + ", " 
            + str(positionsr[i]) + ", "
            + str(positions_refr[i]) + ", "
            + str(velocitiesr[i]) + ", " 
            + str(iqsr[i]) + ", "
            + str(vqsr[i]) + ", " 
            + str(vsupplysr[i])
            + "\n")
f.close()

pygame.quit()

sys.exit()