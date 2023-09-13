import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as ani

I=complex(0,1)

A=0.8       # Amplitude de l'onde
Do=.2       # Largeur initiale de l'enveloppe
ko=10       # Nombre d'onde moyen

# Relation de dispersion (à l'ordre de 2) : 
# w(k) = vP*ko + vG*(k-ko) + vE/2*(k-ko)^2

vP=1.5     # vitesse de phase  (ordre 0)
vG=1       # vitesse de groupe (ordre 1)
vE=0.1     # élargissement     (ordre 2)
wE = vE/Do**2

to=0
xo=-0.3   # point de départ du point bleu (vitesse de phase)

dt=0.005
xmin=-0.5
xmax=1.5

Nbx = 1000

x = np.linspace(xmin, xmax, Nbx)
y = np.linspace(0, 0, Nbx)


def enveloppe (x,t) :
    DDx = Do**2 + I*vE*t
    return A*np.exp(-(x-vG*t)**2/(2*DDx))/np.sqrt(2*np.pi*DDx)

"""
def init(): 
    return 0
"""    
    
def onde (i) :
    t = to + i * dt

    Yenv = enveloppe(x,t)
    Yonde = (Yenv*np.exp(I* 2 * np.pi * (ko*(x-vP*t)))).real
    
    BF.set_ydata (np.abs(Yenv))
    CF.set_data (x,-np.abs(Yenv))
    HF.set_data (x,Yonde.real)
    MBF.set_data (vG*t,np.abs(enveloppe(vG*t,t)))

    a = wE*t/(1+(wE*t)**2)/(2*Do**2)
    b = 2*np.pi*ko - 2*a*vG*t
    c = -2*np.pi*ko*xo - 2*np.pi*ko*vP*t - 0.5*np.arctan(wE*t) + a*(vG*t)**2
    xp = (-b + np.sqrt(b*b-4*a*c))/(2*a)       
    MHF.set_data (xp,np.abs(enveloppe(xp,t)))
 
    affiche_date.set_text (texte_date%(i))


fig = plt.figure ()
ax = fig.add_subplot (111)
BF,  = ax.plot (x, y, lw=2, color='red')
CF,  = ax.plot ([], [], lw=2, color='red')
HF,  = ax.plot ([], [], color='blue')
MBF, = ax.plot ([], [], 'ro', ms=10)
MHF, = ax.plot ([], [], 'bo', ms=7)
ax.set_ylim (-2, 2)
ax.set_xlim (xmin, xmax)

texte_date = 'Date = %.0f'    # affichage du temps
affiche_date = ax.text (0.05, 0.9, '', transform = ax.transAxes)

"""
anim = ani.FuncAnimation (fig, onde, init_func=init, blit=True,\
                                interval=10, repeat=True,frames=300)
"""
anim = ani.FuncAnimation (fig, onde, blit=False,\
                                interval=10, repeat=True,frames=300)

#anim.save ('paquet.mp4', fps = 24)
plt.show()

