import numpy as np 
import matplotlib.pyplot as plt
import iminuit
from iminuit import Minuit
from iminuit.cost import LeastSquares




# [i][j] : i: coincidence number ; j: 0=Lampetia, 1=Nour, 2=Tirgan

def torikomi(file) :
    n=0
    res=[]
    line=file.readline()
    while not line=='' :
        file.readline()
        res.append([0]*3)
        for a in range(3) :
            line=file.readline()
            i = int(line[35]) // 2
            deb = line.find('q1')
            fin = line.find(' ',deb)
            res[n][i]=line[deb+3:fin]
        file.readline()
        line=file.readline()
        n+=1

    file.close()
    res=np.array(res,dtype=int)
    return res

# source=Na22, L1 N3 T4

File = [ '40cm 120deg/','15cm 120deg/','15cm 60deg/' ]
print(File)
ifile= int(input('index?'))

source = open(File[ifile]+'cleaned.txt','r')
bruit  = open(File[ifile]+'cleaned_noise.txt','r')


DATA=torikomi(source)
BRUIT=torikomi(bruit)




# calibration==================

Name=['Lampetia','Nour','Tigran']
Color=['royalblue','darkorange','pink']
Energie = np.array([511,662,1275])


# [i][j] : i: 0=Lampetia, 1=Nour, 2=Tigran ; j: 0=511keV, 1=662keV, 2=1275keV

Channel = np.array([[1.41341e+05,1.81299e+05,3.40647e+05],
                    [1.27389e+05,1.61377e+05,2.96423e+05],
                    [5.40413e+04,7.23864e+04,1.36578e+05]], dtype=float)

Sigma = np.array([  [5.35870e+03,5.77347e+03,7.83329e+03],
                    [4.96903e+03,5.27281e+03,7.45686e+03],
                    [2.33405e+03,2.45340e+03,3.31868e+03]], dtype=float)

Channel_err = np.array([[5.62155e+00,1.02795e+01,1.96977e+01],
                        [4.46016e+00,1.08981e+01,1.37693e+01],
                        [2.06948e+00,5.13651e+00,6.36931e+00]], dtype=float)

Sigma_err = np.array([  [6.38681e+00,9.08245e+00,2.46799e+01],
                        [4.32449e+00,1.19165e+01,1.12525e+01],
                        [2.04232e+00,5.63460e+00,5.72762e+00]], dtype=float)

def linear (x,A,B) :
    return A*x+B


par=[[]]*3

for i in range(3) :
    m=Minuit( LeastSquares(Energie, Channel[i], Sigma[i], linear) , A=100 , B=0 )
    m.migrad()
    a= 1 / m.values['A'] 
    b= -a*m.values['B']
    DATA[:,i]=np.vectorize(linear)(DATA[:,i],a,b)
    if len(BRUIT) :
        BRUIT[:,i]=np.vectorize(linear)(BRUIT[:,i],a,b)
    par[i]=[a,b]

    m=Minuit( LeastSquares(Energie, a*Sigma[i], a*Sigma_err[i], linear) , A=1 , B=0 )
    m.migrad()
    a=m.values['A'] 
    b=m.values['B']
    par[i].append(a)
    par[i].append(b)

par=np.array(par)



def sigma (i,E) :
    return np.vectorize(linear)(E,*par[i][2:])

"""
plt.figure()
for i in range(3) :
    plt.plot(Energie, par[i][0]*Sigma[i] ,'o',color=Color[i])
    plt.plot(Energie, sigma(i,Energie) , 'k-')
plt.show()
"""

#=================================


"""
fig,ax=plt.subplots(3,1)
for i in range(3) :
    ax[i].hist(DATA[:,i],bins=100,range=[0,2000],label=Name[i])
plt.savefig('3gamma.pdf')
"""

#===================================

energie_tot = 1022.
energie_err = np.sum([sigma(i,energie_tot) for i in range(3)])
print(f'energie total = 1022 +/- {energie_err:.0f} keV')

def tri (X):
    res=[]
    for coincidence in X :
        somme = np.sum(coincidence)
        if ( energie_tot - energie_err  <= somme <= energie_tot + energie_err ):
            #for i in range(3):
            #    res.append(coincidence[i])
            res.append(coincidence)
    return np.array(res)


"""
fig,ax=plt.subplots(3,1)
for i in range(3) :
    ax[i].hist(DATA_tri(DATA)[:,i],bins=100,range=[0,2000],label=Name[i])
plt.show()
"""


n_DATA=len(DATA)
n_BRUIT=len(BRUIT)


if input('tri? [y/n]') == 'y' :
    DATA=tri(DATA)
    BRUIT=tri(BRUIT)
    
    print('données : ',len(DATA),'/',n_DATA)
    print('bruit   : ',len(BRUIT),'/',n_BRUIT)





Index=[ [ 0,1,0 ],
        [ 1,2,2 ]]


def gauss (x,a,m,s) :
    return a*np.exp( - (x-m)**2/(2*s**2) )


Energie_coincidence = ([[energie_tot/3]*3]*2+[[273.84407466,273.84407466,474.31185067]])[ifile]

A=[10,60,50]



fig, ax=plt.subplots(3,3)
for i in range(3) :
    i1=Index[0][i]
    i2=Index[1][i]

    ax[0][i].hist2d( DATA[:,i1], DATA[:,i2], range=[[20,1400],[20,1400]], bins=100 , cmap='Greys')
    ax[0][i].set_xlabel(Name[i1]+' (keV)')
    ax[0][i].set_ylabel(Name[i2]+' (keV)')

    ax[0][i].axhline(y=Energie_coincidence[i], c='r' , alpha=.5)   
    ax[0][i].axvline(x=Energie_coincidence[i], c='r' , alpha=.5)

  
    ax[1][i].hist(DATA[:,i],bins=100,color=Color[i],range=[20,1400] ,alpha=1)
    ax[1][i].set_xlabel(Name[i]+' (keV)')

        

    #ax[1][i].axvline(x=energie_coincidence, c='r' , alpha=1)
    X=np.linspace(0,1400,200)
    ax[1][i].plot(X,gauss(X,A[ifile],Energie_coincidence[i],sigma(i,Energie_coincidence[i])),'r',alpha=.5)

    if len(BRUIT) : 
        ax[2][i].hist(BRUIT[:,i],bins=100,color=Color[i],range=[20,1400])
    
    
plt.show()




