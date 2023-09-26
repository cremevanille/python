import numpy as np
import matplotlib.pyplot as plt
from random import random
from time import time
import sys

t0=time()



f=400000*0.9*0.005            #frequence d'emission de 3 gammas 
a=5.6                         #diametre detecteur (cm)
r=float(input('r? (cm)'))     #distance detecteur (cm)

det=int(input('det? [120/60]'))
if det == 120 :
    alphas=[0.,2*np.pi/3,4*np.pi/3]      #angles des detecteurs
else :
    alphas=[0.,5*np.pi/6,7*np.pi/6]
d_alpha = np.arctan( a/(2*r) )

f = f*(d_alpha)**2/2

print((d_alpha)**2/2)


print(f)

spectres=[[],[],[]]

N=int(f*float(input('t? (s)')))


E=1022.

def energies(phi,psi) :
        
    res=[0]*3
            
    res[0] = np.sin(phi + psi) - np.sin(phi)
    res[1] = np.sin(phi + psi)
    res[2] = -np.sin(phi)

    res = res / ( np.sin(psi)*(np.cos(phi) - 1) + np.sin(phi)*(np.cos(psi) - 1) )
    
    res[0] = 1 - res[0]
        
    res = E * res
            
    return res




thetas=np.array([0.]*3)              #angles des gammas


pourcentage=0
ETA=0
n_detection=0

def time_to_string (t) :
    if t > 3600 :
        t = str(round(t//3600))+'h'+str(round(t%3600//60))+'m'+str(round(t%3600%60))+'s'
    elif t > 60 :
        t = str(round(t//60))+'m'+str(round(t%60))+'s'
    else :
        t = str(round(t))+'s'
    return t
    



for n in range(N) :

    if int(n/N*100) > pourcentage :
        pourcentage+=1
        t=time()-t0
        ETA=100*t/pourcentage-t

        bar=str(pourcentage)+'% ['+pourcentage//2*'#'+(50-pourcentage//2)*' '+'] ETA '+time_to_string(ETA)
        print(bar, end='\r')
        sys.stdout.write("\033[K")


    thetas[0]=2*np.pi*random()
    thetas[1]=2*np.pi*random()

        


    x=0
    while x==0 :
        x=random()

    if thetas[0]==thetas[1] :
        
        thetas[2]=thetas[0]+np.pi
        
    
    else :
        
        delta=thetas[1]-thetas[0]

        if 0 <= delta and delta < np.pi : 
            thetas[2]=thetas[0]+np.pi+x*delta
        elif -np.pi < delta and delta < 0 :
            thetas[2]=thetas[1]+np.pi-x*delta
        elif np.pi <= delta and delta < 2*np.pi :
            thetas[2]=thetas[1]+np.pi+x*(2*np.pi-delta)
        else :
            thetas[2]=thetas[0]+np.pi+x*(2*np.pi+delta)
    
    if thetas[2]>=2*np.pi :
        thetas[2]-=2*np.pi

    """
    plt.figure()

    X=np.cos(thetas)
    Y=np.sin(thetas)
 
    for i in range(3) :
        plt.plot([0, X[i]], [0,Y[i]],)
        plt.xlim(-1,1)
        plt.ylim(-1,1)

    plt.show()
    """

    b=np.zeros((3, 3))

    for i in range(len(alphas)) :
        for j in range(len(thetas)) :
            b[i,j] = alphas[i] - d_alpha < thetas[i] and thetas[i] < alphas[i] + d_alpha
    b_detection = bool(np.product(np.array([np.sum(b[i,:]) for i in range(3)])))
    
    if b_detection : 
        # print(np.degrees(thetas))
        
        n_detection+=1


        
        phi=thetas[1] - thetas[0]
        psi=thetas[2] - thetas[1]
        
        #print(energies(phi,psi))
        #print(np.sum(energies(phi,psi)))

        ener=energies(phi,psi)

        for i in range(3) :
            spectres[i].append(ener[i])

bar='100% ['+50*'#'+'] TET '+time_to_string(time()-t0)
print(bar)

print (str(n_detection)+' / '+str(N)+'     '+str(n_detection/N*1000)+'‰')


fig,ax=plt.subplots(3,1)
fig.suptitle('r='+str(r)+' cm // N='+str(N)+' // '+str(n_detection/N*1000)+'‰')
for i in range(3) :
    ax[i].hist(spectres[i],bins=50)




plt.show()
print(energies(5/6*np.pi,np.pi/3))
