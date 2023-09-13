
# -*- coding: utf-8 -*-
"""
@author: Guillaume DELANNOY

Ce document est partagé sous la licence Creative Commons BY-NC-SA : 
Attribution - Pas d’Utilisation Commerciale - Partage dans les Mêmes Conditions.
"""

from __future__ import division

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation

#==============================================
I = complex(0, 1)

periode_animation = 2
images_par_seconde = 50

xmin = -3
xmax = 3

nbx = 500   

def onde(u):
    return np.exp(2*np.pi*I*u)

#==============================================

c = 3e8
gamma = 6e7
epsilon = 9e-12 
frequence = 5e12   

c = 1
epsilon = 1 
frequence = 1

gamma = 1000

#==============================================

w = 2*np.pi * frequence

dt = 1 / (c * periode_animation * images_par_seconde)

Z = (1-I)/np.sqrt(2) * np.sqrt(gamma/(epsilon*w))

rr = (1-Z) / (1+Z)
tt = 2 / (1+Z)

X = np.linspace(xmin, xmax, nbx)

X0=[0,0]
Y0=[-2,2]

Xm = [xmin,xmax]
Ym = np.array([1,1])

fig = plt.figure(figsize=(20, 12))

plt.subplot(211) 
lineV = plt.plot([], [], linewidth=2, color='blue')[0]
plt.plot(X0, Y0, linewidth=3, color='red')
plt.title(u'Champ électrique')
plt.xlim(xmin, xmax)
plt.ylim(-2,2)
plt.grid()

plt.subplot(212) 
lineP = plt.plot([], [], linewidth=2, color='green')[0]
plt.plot(X0, Y0, linewidth=3, color='red')
plt.title(u'Champ magnétique')
plt.xlim(xmin, xmax)
plt.ylim(-2,2)
plt.grid()

def limites():
    plt.plot(Xm, [0,0], linewidth=0.3, color='black')
    plt.plot(Xm, Ym, linewidth=0.3, linestyle='--', color='black')
    plt.plot(Xm,-Ym, linewidth=0.3, linestyle='--', color='black')


def initVP():
    plt.subplot(211) ; limites()
    plt.subplot(212) ; limites() 
    return lineV, lineP,  #il faut les deux sinon erreur au lycée mais pas à la maison !

# mettre un fond de couleur différente pour chaque milieu ?

def animateVP(i): 
    t = i * dt
    Electrique = (onde(c*t-X)+rr*onde(c*t+X)).real*(X<=0) + (tt*onde(c*t-Z*X)).real*(X>=0)
    Magnetique = (onde(c*t-X)-rr*onde(c*t+X)).real*(X<=0) + (Z*tt*onde(c*t-Z*X)).real*(X>=0)
    lineV.set_data(X, Electrique)
    lineP.set_data(X, Magnetique)
    lineP.set_lw(2)
    lineP.set_color('green')
    return lineV, lineP,
 
# pour avoir une image fixe
def statique():
    initVP()     
    t = 0 #51
    to = .0
    Electrique = (onde(c*t-X)+rr*onde(c*t+X)).real*(X<0) + (tt*onde(c*(t-to)-Z*X)).real*(X>0)
    Magnetique = (onde(c*t-X)-rr*onde(c*t+X)).real*(X<0) + (Z*tt*onde(c*(t-to)-Z*X)).real*(X>0)
    plt.subplot(211) ; plt.plot(X,Electrique,'b') ; plt.plot(X,Electrique)
    plt.subplot(212) ; plt.plot(X,Magnetique,'b') ; plt.plot(X,Magnetique)
#    plt.savefig('foo.png', bbox_inches='tight')
    return None

aniVP = animation.FuncAnimation(fig, animateVP, init_func=initVP, frames=periode_animation*images_par_seconde, blit=True, interval=1000//images_par_seconde, repeat=True)

#statique()

plt.show()




