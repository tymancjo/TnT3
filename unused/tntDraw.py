import pygame
import numpy as np
import sys

def init():
    pygame.init()
    global w, h, screen
    w = 800
    h = 800
    screen = pygame.display.set_mode((w, h))

    return screen, w, h
def dot(pos, color, scr, size=5, sx=0, sy=0):
    pygame.draw.circle(scr, color, (int(pos[0]+sx), int(-pos[1]+sy)), size)

def linia(line, color, scr, sx=0, sy=0):

    projection = np.array([(1,0,0),
                            (0,1,0)])

    shiftV = np.array([sx, sy])
    mpl = np.array([1,-1])

    pA = np.matmul(projection, line[0])*mpl + shiftV  
    pB = np.matmul(projection, line[1])*mpl + shiftV

    pA = (int(pA[0]),int(pA[1]))
    pB = (int(pB[0]),int(pB[1]))
    pygame.draw.line(scr, color, pA, pB, 3 )

def draw(bck=(0,0,0)):
    event = pygame.event.poll()
    if event.type == pygame.QUIT:
        sys.exit()
    screen.fill(bck)

def update():
    pygame.display.flip()
    

# while running:
#     draw()
    
#     screen.fill((0, 0, 0))
#     linia(((0,0,0),(100,100,4)),(128,128,128),screen,w/2,h/2)
    
#     update()

    