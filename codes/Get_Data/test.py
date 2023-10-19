import time
import pygame

dt = 0.001

# Inicialização do Pygame
pygame.init()

# Configurações da janela
WIDTH, HEIGHT = 600, 400
screen = pygame.display.set_mode((WIDTH, HEIGHT))
pygame.display.set_caption("Simulação de Temperatura")

# Cores
WHITE = (255, 255, 255)
BLACK = (0, 0, 0)
FONT_COLOR = BLACK

# Fonte
font = pygame.font.Font(None, 36)

# Função para renderizar o texto na tela
def render_text(text, x, y):
    text_surface = font.render(text, True, FONT_COLOR)
    text_rect = text_surface.get_rect()
    text_rect.center = (x, y)
    screen.blit(text_surface, text_rect)


# Mapeamento de teclas para valores de ponto flutuante
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

simulating = True
text = ''

t = time.time()
kp, ki = 0, 0

while simulating:


    screen.fill(WHITE)
    
    events = pygame.event.get()
    for event in events:
        if event.type == pygame.QUIT:
            simulating = False
        elif event.type == pygame.KEYDOWN:
            if event.key == pygame.K_s:  # Pressionar a tecla "s" para finalizar
                simulating = False
            elif event.key == pygame.K_RETURN:
                pid0.Kp = float(textinput.update(events))
                screen.blit(textinput.surface, (10, 10))

            if event.key in key_to_number:
                float_value = key_to_number[event.key]
                print(f"Tecla {pygame.key.name(event.key)} pressionada. Valor float: {float_value}")
                text += float_value

            if event.key == pygame.K_p:
                kp = float(text)
                text = ''
            if event.key == pygame.K_i:
                ki = float(text)
                text = ''
            if event.key == pygame.K_d:
                kd = float(text)
                text = ''



            # Blit its surface onto the screen

    render_text(f"Teste: {text}", WIDTH // 2, HEIGHT // 2)
    render_text(f"Kp: {kp}", WIDTH // 2, HEIGHT // 2 + 30)
    render_text(f"Ki: {ki}", WIDTH // 2, HEIGHT // 2 + 60)

    pygame.display.flip()

    t +=dt
    while(time.time()-t<dt):
        pass


pygame.quit()