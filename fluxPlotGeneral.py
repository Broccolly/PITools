## Import relevant modules

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as anim
import matplotlib as mpl
from tools import *
import f90nml as nml
import os
import scipy
from scipy.signal import hanning, blackman, hamming

plt.rcParams['text.usetex']=True
flip=False
transL=0

#filenames=['./HiRes/Tracking/']
#filenames=['./HiRes/Trb5_A/L00005_000S00_400/']


#filenames=['./HiRes/Tracking/trkBL00050_000S00_500P01/LowS/'] #RPO (a)
#filenames=['./HiRes/Tracking/trkFL00050_000S00_500P01/MidS/'] #RPO (b)
#filenames=['./HiRes/Tracking/trkFL00050_000S00_500P01Cont/HighS/'] #RPO (c)
#filenames=['./HiRes/Tracking/trkFL00050_000S00_500P01Cont/Bif0/'] #RPO (d)
#filenames=['./HiRes/Tracking/trkFL00050_000S00_500P01Cont3/Stable/'] #RPO (e)
#filenames = ['./HiRes/RPO50_0/L00050_000S01_500P85/'] #RPO (f)
#filenames = ['./HiRes/RPO50_0/L00050_000S01_500P83/'] #RPO (g)
#filenames = ['./HiRes/RPO50_0/L00050_000S01_500P81/'] #RPO (h)

#filenames = ['./HiRes/Tracking/SShift50_21/30/']
filenames = ['../PI_RPO/RPO1_478L/POPIOut/']
#filenames = ['../PIModel/BisectOut/']

for j,f in enumerate(filenames):
    states=np.loadtxt(f+'/states.out')
#    plotinfo=nml.read(f+'/PIRun.in')['runparams']
    plotinfo=nml.read(f+'/POPI.in')['popiparams']
    #amps=np.loadtxt(f+'/amps.dat')[-1,:]
    L = plotinfo['L']
#    Nx= plotinfo['N_x']
    Nx= plotinfo['Nx']
    S = plotinfo['S']
    x = np.arange(Nx)*(L/Nx)

    print('S = ',S)
    dt=plotinfo['dt']*plotinfo['spf']
    print('dt = ', dt)
    Nt=states.shape[0]-2
    print(Nt)
    t = np.arange(Nt+1)*dt
    print(L,S)
    tbc=2*np.pi/(L*S)
    print(f, dt, tbc,S)
    #print(states)
    #print(states[0])
    lim=len(t)
    print(Nt)
    uLim=Nt  # Upper and lower limits on timestep number
    lLim=1

    flux=FluxPlot(x[:], t[lLim:uLim], states[lLim:uLim], tbc=tbc, steps=20, title=' ', transL=transL,flip=flip)#, save=f'./Plots/SShiftFlux.png')
    
    fig, ax0=plt.subplots(figsize=(6,2),dpi=200)   
    fig2,ax2=plt.subplots(figsize=(6,2))
    
    y=np.linspace(0,2.*np.pi,100)
    X,Y=np.meshgrid(x,y)
    
    b = 1

    n=states[b,0:Nx]
    E=states[b,Nx:2*Nx]
    phit=states[b,2*Nx:4*Nx:2]+1j*states[b,2*Nx+1:4*Nx:2]
    nt=states[b,4*Nx:6*Nx:2]+1j*states[b,4*Nx+1:6*Nx:2]
    
    # Shift box horizontally (e.g to center structure)
    n=np.roll(n,transL)
    nt=np.roll(nt,transL)
    E=np.roll(E,transL)
    phit=np.roll(phit,transL)
    
    if flip:
        n=-np.flip(n)
        E=-np.flip(E)
        nt=-np.flip(np.conjugate(nt))
        phit=-np.flip(np.conjugate(phit))
    wave=np.exp(1j*Y)
    for x_ in np.arange(np.size(x)):
        wave[:,x_]=wave[:,x_]*nt[x_]+n[x_]
    print(Y)
        
    ax0.plot(x,n,label=r'$\bar{n}$',color='blue')
    ax0.plot(x,E,label=r'$E$',color='red')
    ax0.plot(x,np.abs(nt),label=r'$|\tilde{n}|$'   ,color='blue',linestyle='--')
    ax0.plot(x,np.abs(np.gradient(phit,L/Nx)),label=r'$|\partial_x\tilde{\phi}|$',color='red',linestyle='--')
    ax0.legend()
    ax0.set_title('RPO (e)',pad=-14)
    ax0.set_xlabel(r'$x$',labelpad=-10)
    ax0.set_xlim((x[0],x[-1]))
    fig.tight_layout()

    #ax2.pcolormesh(X,Y,np.real(wave))
    #fig2.tight_layout()
    
    fig2=plt.figure()
    ax2.plot(np.average(flux,axis=1))
    #fig3=plt.figure()
    #plt.plot(np.max(flux,axis=1))
    #plt.savefig('./Plots/SShiftFlux.png')
        #plt.savefig('./Plots/TrbFlux1_5.eps')    

plt.show()
