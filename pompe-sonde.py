import os
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
# plt.rcParams['mathtext.fontset'] = 'cm'
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit


transitions=[1970,1985]



# precision de l'affichage des résultats du fit
nbpix=5 

# temps initial du fit (ps) 
t0s=[5.,7.5,10.,12.5,15.,17.5,20.]

# precision de l'affichage des résultats du fit 
acc=6

# définition de la courbe de fitting
def fit (t,I0,A,T) :
    return I0+A*np.exp(-t/T)

# ordre de grandeur des paramètres attendus
guess=[0.,0.,100.]



plot=int(input('plot? [0/1]\n'))
# plot=0
if plot==1 :
    t0g=float(input('t0? '+str(t0s)+'\n'))
    dexp=int(input('double exp? [0/1]\n'))
tex=int(input('tex? [0/1]\n'))

path='/home/jeremy/Desktop/stage/data/pompe-sonde/'
scan=os.scandir(path)


temps=[]
files=[]

for folder in scan :
    date=folder.name
    scan2=list(os.scandir(path+date))
    for ifile in range(len(scan2)) :

        entry=scan2[ifile].name

        # LECTURE DES FICHIERS

        filepath=path+date+'/'+entry+'/'+entry
        filename=date+'/'+entry

        # detection de la molécule
        if entry[0] == 'W' :
            mol='$\mathrm{W(CO)_6}$'
            # transitions = ...
        else :
            mol='$\mathrm{Fe(CO)_5}$'


        ipol=entry.find('pol')        
        if ipol==-1 :
            pol = '?'
        else :
            ius=entry[ipol:].find('_')
            ib=entry[ipol:].find('b')
            if ib >= 0 : ius=ib
            pol = entry[ipol:][3:ius]


        # lecture des donnés
        delays=np.genfromtxt(filepath+'_delays.csv')
        wavenumbers=np.genfromtxt(filepath+'_wavenumbers.csv')
        signal=np.transpose(np.genfromtxt(filepath+'_signal_sp0_sm0_du0.csv',delimiter=','))

        # on ajoute une entrée pour ce fichier
        temps.append(2*len(t0s)*len(transitions)*['0'])
        files.append(filename)

        
        if plot==1 :

            # PLOT POMP-SONDE 2D 
            plt.figure(filename,figsize=(18,6))
            plt.suptitle(mol+'   '+date[6:]+'/'+date[4:6]+'/'+date[:4]+'   '+entry[-6:-4]+':'+entry[-4:-2]+"'"+entry[-2:]+"''"+'   pol='+pol+'°',size=15)


            plt.subplot(1,len(transitions)+1,1)
            # plt.title('')
            Imax= max(np.abs(signal.min()),signal.max())
            plt.pcolormesh(delays,wavenumbers,signal,cmap='bwr',vmin=-Imax,vmax=Imax)
            plt.xlabel('$t_2~/~\mathrm{ps}$',size=15)
            plt.ylabel('$\omega_3~/~\mathrm{cm^{-1}}$',size=15)
            c=plt.colorbar()
            # c.set_label('intensite, $I~/~\mathrm{u.a.}$',rotation=270)


        # COURBES DE DÉCROISSANCE

        for itps in range(len(t0s)) :
            t0=t0s[itps]

            # indice pour t2=t0
            i0=np.abs(delays-t0).argmin()


            # indice pour t2=1ps
            i1=np.abs(delays-1).argmin()



            for itrans in range(len(transitions)) :
                

                # indice pour w3=transition
                index=np.abs(wavenumbers - transitions[itrans]).argmin()

                # intégration autour de w3=transition
                courbe=np.sum(signal[index-nbpix//2:index+nbpix//2+1],axis=0)/nbpix 

                # calcul du fit apres t0
                par1,cov1=curve_fit(fit,delays[i0:],courbe[i0:],p0=guess)
                err1=np.sqrt(np.diag(cov1))
                courbefit1=np.vectorize(fit)(delays[i0:],*par1)


                if dexp==1 :
                    # calcul du fit avant t0
                    par2,cov2=curve_fit(fit,delays[i1:i0],courbe[i1:i0],p0=guess)
                    err2=np.sqrt(np.diag(cov2))
                    courbefit2=np.vectorize(fit)(delays[i1:i0],*par2)


                # on remplit la case correspondant à cette transition
                j=2*(len(transitions)*itps+itrans)
                temps[ifile][j]=str(par1[2])[:acc]
                temps[ifile][j+1]=str(err1[2])[:acc]


                # PLOT

                if plot==1 and t0g==t0 :

                    ax=plt.subplot(1,len(transitions)+1,itrans+2)

                    # coloration et position text
                    if courbe[i0]>0 :
                        color='r'
                        anchor='upper right'
                        bboxleg=(1,.775)
                        bboxat2=(.638,1)
                    else :
                        color='b'        
                        anchor='lower right'
                        bboxleg=(1,.225)
                        bboxat2=(.638,0)
                    

                    # affichage des paramètres du fit
                    # at=AnchoredText('$I\simeq I_0+A\cdot\exp\left(-\\frac {t}{T}\\right)$\n$I_0='+str(par[0])[:acc]+'\pm'+str(err[0])[:acc]+'~\mathrm{u.a.}$\n$A='+str(par[1])[:acc]+'\pm'+str(err[1])[:acc]+'~\mathrm{u.a.}$\n$T='+str(par[2])[:acc]+'\pm'+str(err[2])[:acc]+'~\mathrm{ps}$\n$t_0='+str(t0)+'~\mathrm{ps}$',frameon=True, loc=anchor)
                    at1=AnchoredText('$I\simeq I_0+A\cdot\exp\left(-\\frac {t-t_0}{T}\\right)$\n$t>t_0='+str(t0)+'~\mathrm{ps}$\n$I_0='+str(par1[0])[:acc]+'~\mathrm{u.a.}$\n$A='+str(par1[1])[:acc]+'~\mathrm{u.a.}$\n$T='+str(par1[2])[:acc]+'~\mathrm{ps}$',frameon=True, loc=anchor)

                    at1.patch.set_boxstyle('round', pad=0.,rounding_size=0.2)

                    ax.add_artist(at1)

                    plt.plot(delays,courbe,color,label='data')
                    plt.plot(delays[i0:],courbefit1,'k',linewidth=.8,label='fit')


                    if dexp==1 :
                        at2=AnchoredText('$I\simeq I_0+A\cdot\exp\left(-\\frac {t}{T}\\right)$\n$1<t<t_0='+str(t0)+'~\mathrm{ps}$\n$I_0='+str(par2[0])[:acc]+'~\mathrm{u.a.}$\n$A='+str(par2[1])[:acc]+'~\mathrm{u.a.}$\n$T='+str(par2[2])[:acc]+'~\mathrm{ps}$',frameon=True, loc=anchor, bbox_to_anchor=bboxat2, bbox_transform=ax.transAxes)

                        at2.patch.set_boxstyle('round', pad=0.,rounding_size=0.2)

                        ax.add_artist(at2)

                        plt.plot(delays[i1:i0],courbefit2,'k',linewidth=.8)

                
                   
                    leg=plt.legend(loc=anchor, bbox_to_anchor=bboxleg)
                    # leg.get_frame().set_edgecolor('k')
                    # leg.get_frame().set_alpha(None)
                    leg.get_frame().set(edgecolor='k',alpha=None)
                    plt.title('$\omega_3='+str(transitions[itrans])+'~\mathrm{cm^{-1}}$',size=15)
                    plt.xlabel('$t_2~/~\mathrm{ps}$',size=15)
                    # plt.ylabel('intensite, $I~/~\mathrm{u.a.}$')
                    plt.axhline(y=0, linewidth=.8, color='k')
                    plt.axvline(x=t0, linewidth=.8, color='k')


        if plot==1 :
            # plt.subplot_tool()
            plt.savefig('output/pompe-sonde/'+date+'_'+entry+'.pdf')
            # plt.show()


pd.DataFrame([[temps[i][2*j]+' +/- '+temps[i][2*j+1] for j in range(len(temps[0])//2)] for i in range(len(temps))],index=files, columns=['w3='+str(x)+' & t0='+str(y) for y in t0s for x in transitions]).to_csv('output/pompe-sonde/temps.csv')

if tex==1 :
    f=open('output/pompe-sonde.tex','w')

    scan=os.scandir('/home/jeremy/Desktop/stage/code/output/pompe-sonde/')
    scan=[x for x in scan if x.name[-4:]=='.pdf']
    for file in scan :
        f.write('\\begin{frame}\n\t\\includegraphics[width=10.5cm]{./pompe-sonde/'+file.name+'}\n\\end{frame}\n')
    f.close()


    os.system('cd output/\npdflatex graph.tex >pompe-sondelog.txt')


