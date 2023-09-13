import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
from scipy.fft import fft,fftfreq
from scipy.interpolate import interp2d
from scipy.optimize import curve_fit



# mauvais calibrage du MCT (cm-1)
probeshift=4

# precision de l'affichage des résultats du fit
acc=5


plot=int(input('plot? [0/1]\n'))
# plot=0


# définition du fit du spectre
def fit (w,I0,A,w0,sig) :
    return I0+A*np.exp(-(w-w0)**2/(2*sig**2))

guess=[0.,50.,2000.,100.]


path='/home/jeremy/Desktop/stage/data/2D-IR/'
scan=os.scandir(path)



for folder in scan :
    date=folder.name
    scan2=os.scandir(path+date)

    # ON NE S'INTERESSE AUX FICHIERS '*_processed.txt'
    scan2=[x for x in scan2 if x.name[-5]=='d']

    for file in scan2 :
        entry=file.name 
        # LECTURE DES FICHIERS

        filepath=path+date+'/'+entry
        filename=date+'/'+entry

        if entry[3] == 'W' :
            mol='$\mathrm{W(CO)_6}$'
        elif entry[3]== 'I' :
            mol='$\mathrm{Fe(CO)_5}$'


        data=np.loadtxt(filepath)



        # SPECTRE DU LASER IR

        # data=dataraw[:,1]
        # # data=np.cos(np.arange(100))
        # t1=np.linspace(-1.0,7.0,len(data))    

        # plt.subplot(1,3,1)
        # plt.title('pyro')
        # plt.xlabel('$t_1~/~\mathrm{ps}$')
        # plt.ylabel('$I$ (u.a.)')
        # plt.plot(t1,data)


        # spectre=np.abs(fft(data-np.mean(data))[:len(data)//2])
        # freq=fftfreq(len(data),8e-12/len(data))[:len(data)//2]/3e8*1e-2

        # plt.subplot(1,3,2)
        # plt.title('TF(pyro)')
        # plt.xlabel('$\omega_1~/~\mathrm{cm^{-1}}$')
        # plt.ylabel('$I$ (u.a.)')
        # plt.ylim(-2,30)
        # plt.plot(freq,spectre)

        # plt.subplot(1,3,3)
        # plt.plot(dataproc[:,1],dataproc[:,2])
        # plt.show()


        # 2D-IR

        # index de la premiere valeur nulle
        i1=data[:,0].argmin()
        probe=data[:i1,0]
        
        pump=data[:,1]
        spectre=data[:,2]
        signal=data[:,3:]

        # x,y=np.meshgrid(probe,pump)
        f=interp2d(probe,pump,signal)

        xmin=1970
        xmax=2030
        ymin=1980
        ymax=2030
        x=np.linspace(xmin,xmax,(xmax-xmin)*10+1)
        y=np.linspace(ymin,ymax,(ymax-ymin)*10+1)
        z=f(x+probeshift,y)
        zmax = max(np.abs(z.min()),z.max())


        # Spectre IR

        

        par,cov=curve_fit(fit,pump,spectre,p0=guess)
        err=np.sqrt(np.diag(cov))
        spectrefit=np.vectorize(fit)(pump,*par)


        if plot==1 :


            fig,(ax1,ax2)=plt.subplots(1,2)
            fig.set_size_inches(12,6)

            plt.suptitle(mol+'   '+date[6:]+'/'+date[4:6]+'/'+date[:4]+'   '+entry[-20:-18]+':'+entry[-18:-16]+"'"+entry[-16:-14]+"''",size=15)

            
            # 2D-IR

            im1=ax1.imshow(z,cmap='bwr',vmin=-zmax, vmax=zmax, extent=[xmin,xmax,ymin,ymax],aspect='equal',origin='lower')
            plt.colorbar(im1,ax=ax1,shrink=.7)

            ax1.contour(x,y,z, colors='k',linewidths=.8,linestyles='-',vmin=z.min(), vmax=z.max(), levels=np.sort([x*y for x in [0.9,0.8,0.6,0.3] for y in [z.min(),z.max()]]))
            
            dmin=max(xmin,ymin)
            dmax=min(xmax,ymax)
            D=np.linspace(dmin,dmax,(dmax-dmin)*10+1)
            ax1.plot(D,D,'k--',linewidth=0.8)


            ax1.set_title('Carte 2D IR')
            ax1.set_xlabel('$\omega_3~/~\mathrm{cm^{-1}}$',size=15)
            ax1.set_ylabel('$\omega_1~/~\mathrm{cm^{-1}}$',size=15)


            # SPECTRE IR

            ax2.set_title('Spectre du laser IR',size=15)
            ax2.set_xlabel('$\omega_1~/~\mathrm{cm^{-1}}$',size=15)
            ax2.plot(pump,spectre,'r',label='data')

            # fit
            # at=AnchoredText('$I\simeq I_0+A\cdot\exp\left(-\\frac {\omega-\omega_0}{2\sigma_\omega^2}\\right)$\n$I_0='+str(par[0])[:acc]+'\pm'+str(err[0])[:acc]+'~\mathrm{u.a.}$\n$A='+str(par[1])[:acc]+'\pm'+str(err[1])[:acc]+'~\mathrm{u.a.}$\n$\omega_0='+str(par[2])[:acc]+'\pm'+str(err[2])[:acc]+'~\mathrm{cm^{-1}}$\n$\sigma_\omega='+str(par[3])[:acc]+'\pm'+str(err[3])[:acc]+'~\mathrm{cm^{-1}}$',frameon=True, loc='upper right')
            at=AnchoredText('$I\simeq I_0+A\cdot\exp\left(-\\frac {\omega-\omega_0}{2\sigma_\omega^2}\\right)$\n$I_0='+str(par[0])[:acc]+'~\mathrm{u.a.}$\n$A='+str(par[1])[:acc]+'~\mathrm{u.a.}$\n$\omega_0='+str(par[2])[:acc]+'~\mathrm{cm^{-1}}$\n$\sigma_\omega='+str(par[3])[:acc]+'~\mathrm{cm^{-1}}$',frameon=True, loc='upper right')
            at.patch.set_boxstyle('round', pad=0.,rounding_size=0.2)
            ax2.add_artist(at)

            ax2.plot(pump,spectrefit,'k',linewidth=.8,label='fit')
            leg=plt.legend(loc='upper right', bbox_to_anchor=(1,.76))
            leg.get_frame().set(edgecolor='k',alpha=None)


            plt.savefig('output/2D-IR/'+date+'_'+entry[:-14]+'.pdf')

    
f=open('output/2D-IR.tex','w')

scan=os.scandir('/home/jeremy/Desktop/stage/code/output/2D-IR/')
scan=[x for x in scan if x.name[-4:]=='.pdf']
for file in scan :
    f.write('\\begin{frame}\n\t\\includegraphics[width=10.5cm]{./2D-IR/'+file.name+'}\n\\end{frame}\n')
f.close()

os.system('cd output\npdflatex graph.tex >2D-IRlog.txt')