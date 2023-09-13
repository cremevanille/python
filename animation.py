#auteur=Jérémy QUÉLIN LECHEVRANTON

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as ani
from math import floor

#plt.rc('text' , usetex=True)

a=-3

def f(n,x) :
    return (n*x)/( 1+n*x**2 )
    #return np.arctan(n*x)
    #return x**a*np.exp(-n*x)*np.sin(x)
    #return n*x**n*np.log(x)

def g(n,x) :
    #return -np.exp(-1)
    #return (a/n)**a*np.exp(-a)
    #return x**a*np.exp(-n*x) 
    return

def frame(i) :
    return (1+10e-3)**i

def l(x) : 
    return 1/x

xmin=-4
xmax=4

ymin=-4
ymax=4

deltat = 8
fps = 60
dt = 1/fps


#==============================================


X=np.linspace(xmin,xmax,1001)


fig = plt.figure()
plt.xlim(xmin,xmax)
plt.ylim(ymin,ymax)
plt.grid(True)




suitefonction1, = plt.plot ([],[], lw=2, color="royalblue" , label='$\\frac{nx}{1+nx^2}$')
suitefonction2, = plt.plot ([],[], lw=2, color="mediumseagreen")
limite = plt.plot (X, l(X) , lw=2, color="crimson" , label='$1/x$')


text = plt.text( xmin*.8 , ymax*.8 , '' , fontsize=12 )
#plt.text(2.5,3 , r'$\frac{1}{x}$' , usetex=True )



def init () :
    suitefonction1.set_data([],[])
    return suitefonction1,

def animate (i) :
    n=frame(i)
    suitefonction1.set_data( X, f(n,X) )
    suitefonction2.set_data( X, g(n,X) )
    text.set_text ('n='+str(floor(10*n)/10))
    return suitefonction1, suitefonction2, text,


anim=ani.FuncAnimation( fig, animate, init_func=init, frames=deltat*fps, interval=1000*dt, blit=True)
plt.legend()
plt.show()

# Writer = ani.writers['ffmpeg']
# writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)



# anim.save('im.mp4', writer=writer)
