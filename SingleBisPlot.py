import numpy as np
import matplotlib.pyplot as plt
import f90nml as nml
import matplotlib.lines as mlines

fname = 'HiRes/Bis5_0/L00005_000S01_200'
fname = '../PIModel/BisectOut/'

amps = np.loadtxt(fname + '/amps.out')
plotInfo = nml.read(fname + '/PIRun.out')
bisInfo = nml.read(fname + '/bisect.in')
states =np.loadtxt(fname + '/states.out')

nframes = np.shape(amps)[1]
nguess = np.shape(amps)[0]

N_t = plotInfo['runparams']['N_t']
N_x = plotInfo['runparams']['N_x']
spf = plotInfo['runparams']['spf']
dt  = plotInfo['runparams']['dt' ]
L   = plotInfo['runparams']['L'  ]
S   = plotInfo['runparams']['S'  ]

turb = bisInfo['bisectparams']['turb']
lam  = bisInfo['bisectparams']['lam']
tmin = bisInfo['bisectparams']['tmin']


phitilde=states[:,2*N_x:4*N_x:2]
phitilde= phitilde +1j*states[:,2*N_x+1:4*N_x:2]
ntilde=states[:,4*N_x:6*N_x:2]
ntilde= ntilde +1j*states[:,4*N_x+1:6*N_x:2]
phiamp = np.absolute(phitilde)
phiamp = np.average(phiamp,axis=1)
namp = np.absolute(ntilde)
namp = np.average(namp,axis=1)



fig = plt.figure(figsize=[8.0,6.0])
t = dt* spf * np.arange(0, nframes)
for g in np.arange(nguess-1):
    ydata = []
    j = 0
    while (amps[g,j] !=0.0 and j < nframes-1):
        ydata.append(amps[g,j])
        j+=1
        if (ydata[-1] > turb):
            col = 'C0'
        else:
            col = 'C1'

    plt.plot(t[0:len(ydata)], ydata, color = col, linewidth = .5, label = 'Attempted Solution')    


ydata = []
j = 0
print(nframes,t[-1])
while (amps[-1,j] !=0.0 and j < nframes - 1):
    ydata.append(amps[-1,j])
    j+=1
    
plt.plot(t[0:len(ydata)], ydata, color = 'C2', linewidth = 2, label='Best Solution')
#plt.plot(t[0:len(ydata)], phiamp[0:-2]/(2*np.pi), color = 'C2', linewidth = .5, label = 'Solution Phase')
#plt.plot(t[0:len(ydata)], namp[0:-2]/(2*np.pi), color = 'C2', linewidth = .5, label = 'Solution Phase')
plt.hlines([lam, turb], 0, t[-1], linestyle = ['-','--'], color='grey')
plt.vlines([tmin],0.5*lam,2*turb, linestyle = ':', color='grey')
#plt.vlines((2*np.pi/(L*S))*np.arange(20),0.5*lam,2*turb, color='purple')

plt.yscale('log')
plt.ylabel(r'$\sqrt{\bar{\tilde{n}^2} +  \bar{\tilde{\phi}^2}}$',fontsize=15)
plt.xlabel('Time',fontsize=15)

turbline  = mlines.Line2D([],[], color = 'C0')
lamline   = mlines.Line2D([],[], color = 'C1')
bestline  = mlines.Line2D([],[], color = 'C2')
lamlim    = mlines.Line2D([],[], linestyle = '-', color = 'grey' )
turblim   = mlines.Line2D([],[], linestyle = '--',  color = 'grey' )
tminlim   = mlines.Line2D([],[], linestyle = ':', color = 'grey' )

#plt.legend([turbline,lamline, bestline, turblim, tminlim],['Eventually\nTurbulent','Eventually\nLaminar', 'Edge State', r'$A_{turb}$', r'$t_{min}$', 'b'], ncol=1,loc='best')
plt.legend([turbline,lamline, bestline, turblim],['Eventually\nTurbulent','Eventually\nLaminar', 'Edge State', r'$A_{turb}$', 'b'], ncol=1,loc='best',fontsize=15)

#plt.ylim((0.5*lam, 1.1*turb))
#plt.ylim((10E-4,1.1*turb))
#plt.ylim((0,0.000001))
#plt.savefig(f'./{fname}/bisect.png')
plt.tight_layout()
plt.show()
